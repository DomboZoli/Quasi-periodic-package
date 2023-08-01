
function [branch,success] = cont_2dtori(varargin)
% perform continuations
%% INPUT
%     varargin:
%       initial cases:
%       1: raw start with initial NR guess: classic situation: u0 (profile+frequencies), par, <cpari>, <pert>
%       2: branch point start: point structure with new parameter set: point, <par>, <cpari>, <pert>
%       3: branch start (from the last point) single point from a branch structure: branch, <par>, <pointN0>, <cpari>, <pert>
%       4: excelling bifurcation stationary solution: u: stationary solution (profile+frequencies+parameter(s)), w: normalized predictor (Dprofile+Dfrequencies+Dparameter(s))
% u, w, par, cpari
%% OUTPUT
%   jc: colocation Jacobian
%% VERSION v13.3
%   v1.0: base code
%   v2.0: maping of numerical derivatives in jac_cond based on cubic interpolation
%   v2.1: allows store convergency details in a point structure by par.num.CONV==1
%   v2.2: modified condition detection
%   v2.3: restricted direction projection can be set (par.num.PR==1 the prediction is in the same direction as it before, bad for FOLD!!!!
%   v2.4: possibility to omit the first NR iteration by using par.num.noINR==1
%   v2.5: try and catch to identify wrong cases
%   v3.0: new jac_cond based on interpolation and patches
%   v3.1: jac_cond: minimum of the minimums of the feasable solutions
%   v3.2: copy of 3.1
%   v4.0: check condition function for multi dof continuation
%   v4.0_f1f2: modified for two freed frequencies, for turning model
%   v4.1: modified initial start extended with perturbation values
%   v4.3: bug in branch initialization has been solved
%   v4.4: analytical parameter derivatives in jac_cpari
%   v4.5: analytical derivatives in jac_cond
%   v4.6: milling case with flyover
%   v5.0: fixing one frequency by par.num.fixom=0, 1 or 2 by truncation on J
%   v5.1: convergence success; storing conv debuged
%   v5.2: plot only real not complex
%   v6.0: introduction of sys_rhs and sys_deri
%   v6.1: selectable scheme, par.num.scheme, 0: equidistant, 1: Chebishev, 2: Legendre
%   v6.2: continuation problems with fixom corrected
%   v7.0: continuation parameter boundaries
%   v7.1: bug at boundary hit
%   v7.2: jac_cond is debugged
%   v8.0: direct excel of a branch from the assumed bifurcation point by using two initial vectors
%   v9.0: idealistic case distinguished from flyover, par.system.NFO==0: ideal, par.system.NFO>0 flyover
%   v9.1: different predictor color in cont_2dtori heading
%   v9.2: dim and dof in the par.system list, def_taus introduced
%   v9.3: debugged bounds
%   v9.4: debugged 2
%   v9.5: debugged 3 delayed continuation
%   v10.0: condition function outsourced
%   v11.0: jac_cpari outsourced to jac_cpari_ex by using numpars
%   v11.1: norm control by checking the solution as well
%   v11.2: representation point shape functions and their derivatives, for condition definition
%   v11.3: plot first the 0th iteration step
%   v12.0: numerical derivation of the condition function
%   v12.1: initialize with different resolution and interpolate profile
%   v12.2: possibiliy to save batch in each step (for anyway slow calc)
%   v12.3: possible to change between phase condition (pc) types: par.num.pctype==0: Poincaré pc, par.num.pctype==1: least square type
%   v12.4: keeping previous solution for phase condition
%   v13.0: relativelly cleared up code
%   v13.1: parameters for added conditions sys_cond condpar is added to numpars defined in par.system.condpar
%   v13.2: tic;toc has added to the list 
%   v13.3: debugged
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% CONSTANTS
fig1=figure(1);
fig2=figure(2);
%plot 1par bifdiagram
PNP=1;
%plot parameter
PPAR=1;
%allowed printed lines
NPL=5;
%position tolerance
uTOL=1e-8;
%predictor color
prec='m.-';
%numerical condition determination
NUMCOND=false;
%stepsave=[] (no) or name stepsavename
stepsavename=[];
%% backup branch
global bkpbranch;
bkpbranch=[];
%% PARAMETERS
%initialy no intial branch
brinit=false;
%found initial case
ficase=false;
iwi=[];
% case study about the different starting options with varargin
%% INITIALIZATION case 1: simple start from a point
if ~isstruct(varargin{1}) && length(varargin)>=2 && isstruct(varargin{2})
    %branch initial
    brinit=false;
    %classic situation: u0 (profile+frequencies), par, <cpari>, <pert>
    u0j=varargin{1};
    %set perturbation
    if length(varargin)>2 && isnumeric(varargin{3})
        pert=varargin(3);
    else
        pert=1;
    end
    cNom1=varargin{2}.num.NC1*varargin{2}.num.p+1;
    cNom2=varargin{2}.num.NC2*varargin{2}.num.p+1;
    %proceed perturbation
    if ~(~isempty(pert) && pert~=1)
        pert=1;
    end
    %set dimension
    n=round((length(u0j)-2)/cNom1/cNom2);
    %faulty initial condition data structure
    ficase=true;
    par=varargin{2};
    if length(varargin)>2
        cpari=varargin{3};
        if length(cpari)>=1
            for kl=1:length(cpari)
                u0j(cNom1*cNom2*n+2+kl)=par.system.(par.parlist{cpari(kl)});
            end
        end
    end
    %continuation step counter
    lcont=0;
end
%% INITIALIZATION case 2
if isstruct(varargin{1}) && isfield(varargin{1},'profile')
    %point structure with new parameter set: point, <par>, <cpari>, <pert>
    cNom1=varargin{1}.par.num.NC1*varargin{1}.par.num.p+1;
    cNom2=varargin{1}.par.num.NC2*varargin{1}.par.num.p+1;
    %set perturbation
    if length(varargin)>3 && isnumeric(varargin{4})
        pert=varargin(4);
    else
        pert=1;
    end
        n=size(varargin{1}.profile,3);
    %proceed perturbation
    if ~(~isempty(pert) && pert~=1)
        pert=1;
    end
    
    if length(varargin)>1 && isfield(varargin{2},'num')
        par=varargin{2};
    else
        %varargin{2}=[]!
        par=varargin{1}.par;
    end
    %stepsavename
    if isfield(par.num,'stepsavename')
        stepsavename=par.num.stepsavename;
    end
    %TODO: here different kind of perturbation must be set
    if ~(cNom1==par.num.NC1*par.num.p+1 && cNom2==par.num.NC2*par.num.p+1)
        nNom1=par.num.NC1*par.num.p+1;
        nNom2=par.num.NC2*par.num.p+1;
                m1=linspace(0,2*pi,cNom1);
        m2=linspace(0,2*pi,cNom2);
        nm1=linspace(0,2*pi,nNom1);
        nm2=linspace(0,2*pi,nNom2);
        %interpolate profile
        nprofile=zeros(nNom1,nNom2,n);
        for l_int=1:n
            nprofile(:,:,l_int)=interp2(m1,m2,varargin{1}.profile(:,:,l_int)',nm1,nm2.','cubic')';
        end
        varargin{1}.profile=nprofile;
%         varargin{1}.mesh.m1=nm1;
%         varargin{1}.mesh.m2=nm2;
        varargin{1}.par.num.NC1=par.num.NC1;
        varargin{1}.par.num.NC2=par.num.NC2;
        varargin{1}.par.num.p=par.num.p;
        cNom1=nNom1;
        cNom2=nNom2;
    end
    u0j=kl2j(pert*varargin{1}.profile,cNom1,cNom2);
    
    u0j(end+1)=varargin{1}.omega1;
    u0j(end+1)=varargin{1}.omega2;
    ficase=true;

    if length(varargin)>1
        if isnumeric(varargin{2})
            cpari=varargin{2};
        else
            cpari=varargin{3};
        end
        if length(cpari)>=1
            for kl=1:length(cpari)
                u0j(cNom1*cNom2*n+2+kl)=par.system.(par.parlist{cpari(kl)});
            end
        end
    end
    %branch initial, no NR iteration is required
    brinit=false;
    %continuation step counter
    lcont=0;
end
%% INITIALIZATION case 3:
if isstruct(varargin{1}) && isfield(varargin{1},'points')
    %single point from a branch structure: branch, <par>, <pointN0>, <cpari>, <pert>
    if length(varargin)>1
        if length(varargin)>2
            pN0=varargin{3};
            lcont=0;
        else
            pN0=length(varargin{1}.points);
            lcont=pN0;
        end
        if isfield(varargin{2},'num')
            par=varargin{2};
        else
            par=varargin{1}.points(pN0).par;
        end
    else
        pN0=length(varargin{1}.points);
        par=varargin{1}.points(pN0).par;
        lcont=pN0;
    end
    %set perturbation
    if length(varargin)>4 && isnumeric(varargin{5})
        pert=varargin(5);
    else
        pert=1;
    end
    cNom1=varargin{1}.points(pN0).par.num.NC1*varargin{1}.points(pN0).par.num.p+1;
    cNom2=varargin{1}.points(pN0).par.num.NC2*varargin{1}.points(pN0).par.num.p+1;
    if ~isempty(pert) && pert~=1
        %TODO: here different kind of perturbation must be set
        u0j=kl2j(pert*varargin{1}.points(pN0).profile,cNom1,cNom2);
    else
        u0j=kl2j(varargin{1}.points(pN0).profile,cNom1,cNom2);
    end
    n=size(varargin{1}.points(pN0).profile,3);
    u0j(end+1)=varargin{1}.points(pN0).omega1;
    u0j(end+1)=varargin{1}.points(pN0).omega2;
    branch=varargin{1};
    ficase=true;
    if length(varargin)>3
        cpari=varargin{4};
        if length(cpari)>=1
            for kl=1:length(cpari)
                u0j(cNom1*cNom2*n+2+kl)=par.system.(par.parlist{cpari(kl)});
            end
        end
        %branch initial
        brinit=false;
    else
        cpari=varargin{1}.contpar;
        if length(cpari)>=1
            for kl=1:length(cpari)
                u0j(cNom1*cNom2*n+2+kl)=par.system.(par.parlist{cpari(kl)});
            end
        end
        %branch initial
        brinit=true;
    end
end
%% INITIALIZATION case 4:
if ~isstruct(varargin{1}) && length(varargin)>3 && ~isstruct(varargin{2}) && isstruct(varargin{3}) && ~isstruct(varargin{4})
    %branch initial
    brinit=false;
    %classic situation: u0 (profile+frequencies), par, <cpari>, <pert>
    u0j=varargin{1};
    %initial of prediction
    iwi=varargin{2};
    u0j=u0j+0.1*iwi/norm(iwi);
    %cpari
    cpari=varargin{4};
    %perturbation, dummy parameter
    pert=1;
    %dimensions
    cNom1=varargin{3}.num.NC1*varargin{3}.num.p+1;
    cNom2=varargin{3}.num.NC2*varargin{3}.num.p+1;
    %dimension of the system
    n=round((length(u0j)-2-length(cpari))/cNom1/cNom2);
    %continuation step counter
    lcont=0;
    if length(iwi)==cNom1*cNom2*n+2+length(cpari)
        %parameters
        par=varargin{3};
        %continuation parameters
        cpari=varargin{4};
        if length(cpari)>=1
            for kl=1:length(cpari)
                u0j(cNom1*cNom2*n+2+kl)=par.system.(par.parlist{cpari(kl)});
            end
        end
        %normalization
%         iwi=iwi/norm(iwi);
    else
        ficase=false;
    end
    ficase=true;
end
par.num.cpari=cpari;
%check initial data
if ~ficase
    error('Faulty initial parameters or function argument structure!');
end
%no NR initial iteration
if par.num.noINR
    brinit=true;
end
%finite differences
switch par.num.dp
    case 2
        etas=[-1/2 0 +1/2];
    case 4
        etas=[1/12 -8/12 0 +8/12 -1/12];
    case 6
        etas=[-1/60 9/60 -45/60  0 +45/60 -9/60 1/60];
    case 8
        etas=[3/840 -32/840  168/840  -672/840 0 672/840 -168/840 32/840 -3/840];
end
%allowed NRI
aNRI=par.num.maxNRI;
% collocation intervals
NC1=par.num.NC1;
NC2=par.num.NC2;
%set default for phase condition
if ~isfield(par.num,'pctype');par.num.pctype=1;end
%polynomial order
P=par.num.p;
% Nom=10;
Nom1=NC1*P+1;
% discretization parameter for the vibration frequency
Nom2=NC2*P+1;
%calculate the delays
[~,nptk]=par.system.sys_tau([],par);
%delay parameter index
dpik=find(nptk==cpari,1);
%the parameter is a delay? if yes G must be stored for cpari
knpisdelayk=~isempty(dpik) || numel(cpari)>1;
%coefficients
cpps=zeros(7,7);
switch par.num.scheme
    case 0
        %equidistant
        cpps(1,1)=0.5;
        cpps(2,1:2)=[1/3 2/3];
        cpps(3,1:3)=[1/4 1/2 3/4];
        cpps(4,1:4)=[1/5 2/5 3/5 4/5];
        cpps(5,1:5)=[1/6 1/3 1/2 2/3 5/6];
        cpps(6,1:6)=[1/7 2/7 3/7 4/7 5/7 6/7];
        cpps(7,1:7)=[1/8 1/4 3/8 1/2 5/8 3/4 7/8];
    case 1
        % Chebyshev constants
        cpps(1,1)=0.5;
        cpps(2,1:2)=[0.146447, 0.853553];
        cpps(3,1:3)=[0.0669873, 0.5, 0.933013];
        cpps(4,1:4)=[0.0380602, 0.308658, 0.691342, 0.96194];
        cpps(5,1:5)=[0.0244717, 0.206107, 0.5, 0.793893, 0.975528];
        cpps(6,1:6)=[0.0170371, 0.146447, 0.37059, 0.62941, 0.853553, 0.982963];
        cpps(7,1:7)=[0.012536, 0.109084, 0.283058, 0.5, 0.716942, 0.890916, 0.987464];
    case 2
        %Legendre
        cpps(1,1)=0.5;
        cpps(2,1:2)=[0.211325, 0.788675];
        cpps(3,1:3)=[0.112702, 0.5, 0.887298];
        cpps(4,1:4)=[0.0694318, 0.330009, 0.669991, 0.930568];
        cpps(5,1:5)=[0.0469101, 0.230765, 0.5, 0.769235, 0.95309];
        cpps(6,1:6)=[0.0337652, 0.169395, 0.38069, 0.61931, 0.830605, 0.966235];
        cpps(7,1:7)=[0.025446, 0.129234, 0.297077, 0.5, 0.702923, 0.870766, 0.974554];
    case 3
        %Lobatto
        cpps(1,1)=0.5;
        cpps(2,1:2)=[1/3 2/3];
        cpps(3,1:3)=[0, 0.5, 1];
        cpps(4,1:4)=[0, 0.276393202250021,  0.723606797749979, 1];
        cpps(5,1:5)=[0,   0.172673164646011, 0.5,   0.827326835353989, 1];
        cpps(6,1:6)=[0,0.117472338035268, 0.357384241759677, 0.642615758240322, 0.882527661964732, 1];
        cpps(7,1:7)=[0,  0.0455369370651239, 0.280263200471634, 0.5,0.719736799528366,       0.954463062934876, 1];
end
%select base order constants
cp=cpps(P,1:P);
%representation points in an interval
tp=linspace(0,1,P+1);
%all shape functions in all collocation points
cPpq=zeros(P,P,P+1,P+1);
tPpq=zeros(P+1,P+1,P+1,P+1);
cd1Ppq=zeros(P,P,P+1,P+1);
cd2Ppq=zeros(P,P,P+1,P+1);
td1Ppq=zeros(P+1,P+1,P+1,P+1);
td2Ppq=zeros(P+1,P+1,P+1,P+1);
clag=zeros(P+1,P);
cdlag=zeros(P+1,P);
tlag=zeros(P+1,P+1);
tdlag=zeros(P+1,P+1);
for kP=0:P
    clag(kP+1,:)=lag(cp,kP);
    cdlag(kP+1,:)=dlag(cp,kP);
    tlag(kP+1,:)=lag(tp,kP);
    tdlag(kP+1,:)=dlag(tp,kP);
end
% collocation interval length
dth1=2*pi/NC1;
dth2=2*pi/NC2;
%interval points
th1i=((1:NC1+1)-1)*dth1;
th2j=((1:NC2+1)-1)*dth2;
%collocation points
c1ip=reshape(repmat(th1i(1:end-1),[P,1])+repmat(cp'*dth1,[1 NC1]),[NC1*P,1])';
c2jp=reshape(repmat(th2j(1:end-1),[P,1])+repmat(cp'*dth2,[1 NC2]),[NC2*P,1])';
for pp=1:P+1
    for qq=1:P+1
        %on collocation points
        cPpq(:,:,pp,qq)=clag(pp,:)'*clag(qq,:);
        cd1Ppq(:,:,pp,qq)=cdlag(pp,:)'*clag(qq,:)/dth1;
        cd2Ppq(:,:,pp,qq)=clag(pp,:)'*cdlag(qq,:)/dth2;
        %on representation points
        tPpq(:,:,pp,qq)=tlag(pp,:)'*tlag(qq,:);
        td1Ppq(:,:,pp,qq)=tdlag(pp,:)'*tlag(qq,:)/dth1;
        td2Ppq(:,:,pp,qq)=tlag(pp,:)'*tdlag(qq,:)/dth2;
    end
end
%define numpars
numpars.Nom1=Nom1;
numpars.Nom2=Nom2;
numpars.P=P;
numpars.n=n;
numpars.cPpq=cPpq;
numpars.cd1Ppq=cd1Ppq;
numpars.cd2Ppq=cd2Ppq;
numpars.tPpq=tPpq;
numpars.td1Ppq=td1Ppq;
numpars.td2Ppq=td2Ppq;
numpars.c1ip=c1ip;
numpars.c2jp=c2jp;
numpars.NC1=NC1;
numpars.NC2=NC2;
%define sys_cond related specific parameters
if isfield(par.system,'condpar')
    numpars.condpar=par.system.condpar;
end

%number of parameters to be continued
if exist('cpari','var')
    if ~brinit
        %to have initial NR iteration
        np=length(cpari)-1;
    else
        %jump to continuation
        np=length(cpari);
    end
else
    np=0;
end
%representation steps
h1=2*pi/(NC1*P);
h2=2*pi/(NC2*P);
%representation meshpoints
mth1=0:h1:(NC1*P)*h1;
mth2=0:h2:(NC2*P)*h2;

%check if interpolation is needed in the initial point
if cNom1~=Nom1 || cNom2~=Nom2
    %interpolation of initial point
    %rest of the parameters stored in v0
    pars0=u0j(cNom1*cNom2*n+1:end);
    %2D representation
    u0kl=j2kl(u0j(1:cNom1*cNom2*n),cNom1,cNom2);
    %size of extension
    eN=3;
    u0kle=zeros(cNom1+2*eN,cNom2+2*eN,size(u0kl,3));
    u0kle(eN+1:end-eN,1:eN,:)=u0kl(:,(end-(eN-1)):end,:);
    u0kle(eN+1:end-eN,(end-(eN-1)):end,:)=u0kl(:,1:eN,:);
    u0kle(eN+1:end-eN,eN+1:end-eN,:)=u0kl;
    u0kle(1:eN,:,:)=u0kle((end-(2*eN-1)):end-eN,:,:);
    u0kle((end-(eN-1)):end,:,:)=u0kle(eN+1:2*eN,:,:);
    %old phase span
    phisi1=linspace(-2*pi/cNom1*eN,2*pi+(eN-1)*2*pi/cNom1,cNom1+2*eN);
    phisi2=linspace(-2*pi/cNom2*eN,2*pi+(eN-1)*2*pi/cNom2,cNom2+2*eN);
    phis1=linspace(0,2*pi,Nom1);
    phis2=linspace(0,2*pi,Nom2);
    %interpolation
    ipukl=zeros(Nom1,Nom2,n);
    for l=1:n
        ipukl(:,:,l)=interp2(phisi1,phisi2,u0kle(:,:,l)',phis1,phis2','spline')';
    end
    u0j=kl2j(ipukl);%+kl2j(normrnd(uTOL,uTOL/2*sqrt(Nom1),[Nom1,Nom2,n]));
    u0j(Nom1*Nom2*n+1:Nom1*Nom2*n+length(pars0))=pars0;
    %NR is required
    brinit=false;
end
%% MAIN CODE
%modal Dof
mn=par.system.dof;
%dimension of the system + 2 dimensional oscillatory system
n=2*mn;
%order of the cental difference
m=(length(etas)-1);
mp=round(m/2);
%bifdiagram
plb=zeros(1+1+np,0);
%printed lines
PL=0;
%break downs
NBD=0;
%overshoots
NOS=0;
%healthy continuation
NHC=0;
%steptype
steptype='in';
%original step length befor break down
dsBD=par.num.ds;

%maximum NR iterations
maxNR=par.num.maxNR;
conver=[];
wpi=[];
upj=zeros(size(u0j));
tic;
%time for one continued point
pcsten=[]; newstep=false;
%% continuation
while (lcont<=par.num.Nc && NOS<=par.num.maxOVS)
    pcstst=tic;
    %start the iteration
    k_nit=1;
    %kill increase norm condition
    kNRI=0;
    %set norm big
    dnorm=1000;
    %% Newton Raphson iteration (correction)
    %set initial condition
    pjp=u0j;
    if ~brinit
        %par.num.conv
        if par.num.CONV==1
            conver=[];
        end
        pitsten=[];
        while (((k_nit<=par.num.maxNR && kNRI<=par.num.maxNRI && lcont==0) || ...
                (k_nit<=par.num.maxNRc && kNRI<=par.num.maxNRIc && lcont>0)) && ...
                dnorm>=par.num.absTOL && sum(isnan(pjp))==0 && sum(isinf(pjp))==0)
            %timer for iteration start
            pitstst=tic;
            %keepint previous solution
            if k_nit~=1;pj0=pj;else pj0=u0j;end
%             if k_nit==1;pj0=u0j;end%else pj0=u0j;end
            %replace iterative solution
            pj=pjp;
            %delays
            taus=par.system.sys_tau([],par);
            %print solution
            if par.num.disp(1) && PL==0
                if np==1
                    disp(['#P' 9 'TR(s)' 9 '#NR' 9 '#NRI' 9 '#NBD' 9 '#NOS' 9 par.parlist{cpari(1)} 9 9 9 'freq_1 (Hz)' 9 'freq_2 (Hz)' 9 9 '||u||' 9 9 '||ru||' 9 9 '||G||']);
                elseif np==2
                    disp(['#P' 9 'TR(s)' 9 '#NR' 9 '#NRI' 9 '#NBD' 9 '#NOS' 9 par.parlist{cpari(1)} 9 9 9 par.parlist{cpari(2)} 9 9 9 'freq_1 (Hz)' 9 'freq_2 (Hz)' 9 9 '||u||' 9 9 '||ru||' 9 9 '||G||']);
                else
                    disp(['#P' 9 'TR(s)' 9 '#NR' 9 '#NRI' 9 '#NBD' 9 '#NOS' 9 9 9 9 9 9 9 'freq_1 (Hz)' 9 'freq_2 (Hz)' 9 9 '||u||' 9 9 '||ru||' 9 9 '||G||']);
                end
                disp('------------------------------------------------------------')
            end
            %determine iteration time
            if ~isempty(pitsten);Ditst=pitsten;else;Ditst=[];end
            %display
            if par.num.disp(1)==1 && par.num.disp(2)==1 && k_nit==1 && lcont==0
                if np==1
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9 num2str(par.system.(par.parlist{cpari(1)}),'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n')]);
                elseif np==2
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9 num2str(par.system.(par.parlist{cpari(1)}),'%10.5e\n') 9 num2str(par.system.(par.parlist{cpari(2)}),'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n')]);
                else
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9  9  9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n') 9 9 num2str(nan,'%10.5e\n')]);
                end
                %                 minths
                PL=PL+1;
            end
            if ~isempty(pitsten);pitsten=[];end
            %building the Jacobian for the torus and one parameter continuation
            J=spalloc(Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+np,ceil(0.1*(Nom1*Nom2*n+2)^2));
            %set solution to be compared
            par.num.uj0=pj0(1:Nom1*Nom2*n,1);
            %define sparse array for Jacobians
            if ~knpisdelayk
                J(1:Nom1*Nom2*n+2,1:Nom1*Nom2*n+2)=jac_collocation_2dtori(pj,par);
            else
                [J(1:Nom1*Nom2*n+2,1:Nom1*Nom2*n+2),outtau,outp]=jac_collocation_2dtori(pj,par);
                out{1}=outtau;
                out{2}=outp;
            end
            %include continuation conditions (only if this is not the first correction)
            if lcont>0
                %derivatives w.r.t. parameters
                if ~knpisdelayk
                    numpars.np=np;
                    J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+3:Nom1*Nom2*n+2+np)=jac_cpari_ex(pj,par,cpari,[],numpars);
                else
                    numpars.np=np;
                    J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+3:Nom1*Nom2*n+2+np)=jac_cpari_ex(pj,par,cpari,out,numpars);
                end
                %set condition
                if np==2
                    numpars.np=np;
                    if ~NUMCOND
                        J(Nom1*Nom2*n+2+1,1:Nom1*Nom2*n+2+np)=par.system.sys_cond(pj,par,cpari,numpars,out);
                    else
                        Ju=zeros(1,Nom1*Nom2*n+2+np);
                         Jk=par.system.sys_cond(pj,par,cpari(1),numpars,out);
                        parfor k=1:Nom1*Nom2*n
                            %parfor
                            pjpp=pj;
                            pjpp(k)=pjpp(k)+uTOL;
                            parp=par;
                            if k>Nom1*Nom2*n+2
                                parp.system.(par.parlist{cpari(k-(Nom1*Nom2*n+2))})=pjpp(k);
                            end
                            condp=par.system.sys_cond(pjpp,parp,[],numpars,out);
                            pjmm=pj;
                            pjmm(k)=pjmm(k)-uTOL;
                            parm=par;
                            if k>Nom1*Nom2*n+2
                                parm.system.(par.parlist{cpari(k-(Nom1*Nom2*n+2))})=pjmm(k);
                            end
                            condm=par.system.sys_cond(pjmm,parm,[],numpars,out);
                            Ju(1,k)=(condp-condm)/(2*uTOL);
                        end
                        Ju(1,Nom1*Nom2*n+1:Nom1*Nom2*n+2+np)=Jk(1,Nom1*Nom2*n+1:Nom1*Nom2*n+2+np);
                        J(Nom1*Nom2*n+2+1,1:Nom1*Nom2*n+2+np)=Ju;
                    end
                end
                if ~instepEB
                    %set distance condition for pseudo arclength (upj: previous solution in the branch)
                    J(Nom1*Nom2*n+2+np,1:Nom1*Nom2*n+2+np)=2*(pj.'-upj.');
                else
                    J(Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+np)=1;
                end
                %calculate Newton-Raphson step
                if par.num.fixom==0
                    itrunc=1:Nom1*Nom2*n+2+np;
                else
                    pjp=pj;
                    %remove one of the phase condition
                    itrunc=[1:Nom1*Nom2*n+par.num.fixom-1 Nom1*Nom2*n+par.num.fixom+1:Nom1*Nom2*n+2+np];
                end
                if ~instepEB
                    numpars.np=np;
                    fi=fun_invariance_cont(pj,par,cpari,upj,0,numpars);
                    pjp(itrunc,1)=pj(itrunc,1)-J(itrunc,itrunc)\fi;
                else
                    numpars.np=np;
                    fi=fun_invariance_cont(pj,par,cpari,upj,1,numpars);
                    pjp(itrunc,1)=pj(itrunc,1)-J(itrunc,itrunc)\fi;
                end
            else
                %only Newton-Raphson of the profile
                if np==0
                    %simple one parameter continuation start
                    if par.num.fixom==0
                        itrunc=1:Nom1*Nom2*n+2;
                    elseif par.num.fixom==1
                        pjp=pj;
                        itrunc=[1:Nom1*Nom2*n Nom1*Nom2*n+2];
                    else
                        pjp=pj;
                        itrunc=1:Nom1*Nom2*n+1;
                    end
                    numpars.np=np;
                    fi=fun_invariance_cont(pj(1:Nom1*Nom2*n+2,1),par,[],[],[],numpars);
                    pjp(itrunc,1)=pj(itrunc,1)-J(itrunc,itrunc)\fi;
                else
                    %fixed omega indeces
                    if par.num.fixom==0
                        itrunc=1:Nom1*Nom2*n+2+np;
                    else
                        %remove one of the phase condition
                        itrunc=[1:Nom1*Nom2*n+par.num.fixom-1 Nom1*Nom2*n+par.num.fixom+1:Nom1*Nom2*n+2+np];
                    end
                    %extension by the first parameter
                    if ~knpisdelayk
                        numpars.np=np;
                        J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+np)=jac_cpari_ex(pj,par,cpari(1),[],numpars);
                    else
                        numpars.np=np;
                        J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+np)=jac_cpari_ex(pj,par,cpari(1),out,numpars);
                    end
                    %extension by condition part
                    numpars.np=np;
                    if ~NUMCOND
                        J(Nom1*Nom2*n+2+1,1:Nom1*Nom2*n+2+np)=par.system.sys_cond(pj,par,cpari(1),numpars,out);
                    else
                        Ju=zeros(1,Nom1*Nom2*n+2+np);
                        Jk=par.system.sys_cond(pj,par,cpari(1),numpars,out);
                        parfor k=1:Nom1*Nom2*n
                            %parfor
                            pjpp=pj;
                            pjpp(k)=pjpp(k)+uTOL;
                            parp=par;
                            if k>Nom1*Nom2*n+2
                                parp.system.(par.parlist{cpari(k-(Nom1*Nom2*n+2))})=pjpp(k);
                            end
                            condp=par.system.sys_cond(pjpp,parp,[],numpars,out);
                            pjmm=pj;
                            pjmm(k)=pjmm(k)-uTOL;
                            parm=par;
                            if k>Nom1*Nom2*n+2
                                parm.system.(par.parlist{cpari(k-(Nom1*Nom2*n+2))})=pjmm(k);
                            end
                            condm=par.system.sys_cond(pjmm,parm,[],numpars,out);
                            Ju(1,k)=(condp-condm)/(2*uTOL);
                        end
                        Ju(1,Nom1*Nom2*n+1:Nom1*Nom2*n+2+np)=Jk(1,Nom1*Nom2*n+1:Nom1*Nom2*n+2+np);
                        J(Nom1*Nom2*n+2+1,1:Nom1*Nom2*n+2+np)=Ju;
                    end
                    %two parameter continuation start
                    if par.num.fixom==0
                        numpars.np=np;
                        fi=fun_invariance_cont(pj(1:Nom1*Nom2*n+2+np,1),par,cpari(1),[],1,numpars);
                        pjp=pj(1:Nom1*Nom2*n+2+1)-J\fi;
                    else
                        %trunction with fixed frequencies
                        numpars.np=np;
                        fi=fun_invariance_cont(pj(1:Nom1*Nom2*n+2+np,1),par,cpari(1),[],1,numpars);
                        pjp(itrunc,1)=pj(itrunc,1)-J(itrunc,itrunc)\fi;
                    end
                end
            end
            %pkl
            pkl=j2kl(pj(1:Nom1*Nom2*n));
            pklp=j2kl(pjp(1:Nom1*Nom2*n));
            %reset parameters calculated by NR iteration
            if length(pjp)>Nom1*Nom2*n+2 && np>0
                par.system.(par.parlist{cpari(1)})=pjp(Nom1*Nom2*n+3);
            end
            if length(pjp)>Nom1*Nom2*n+3 && np>1
                par.system.(par.parlist{cpari(2)})=pjp(Nom1*Nom2*n+4);
            end
            if par.num.plot(1)
                %profile
                if par.num.plot(2)
                    set(0,'CurrentFigure',fig1);
                    hold off;surf(real(pklp(:,:,1)),'edgecolor','none','facecolor','interp');
                    xlabel('\theta_2');
                    ylabel('\theta_1');
                    title(int2str(k_nit));
                    drawnow;
                end
                %infinite norm plot for converging points
                if par.num.plot(3)
                    set(0,'CurrentFigure',fig2);
                    hold on;
                    if length(cpari)==1 && PNP==1
                        plot(par.system.(par.parlist{cpari(PPAR)}),max(max(abs(pklp(:,:,1)),[],1),[],2),'bx');
                    elseif np==0
                        plot(max(max(abs(pklp(:,:,1)),[],1),[],2),max(max(abs(pklp(:,:,2)),[],1),[],2),'bx');
                    else
                        plot(par.system.(par.parlist{cpari(2)}),par.system.(par.parlist{cpari(1)}),'bx');
                    end
                    drawnow;
                end
            end
            %calculate difference norm for the position coordinate
            %             dnormn=norm(reshape(pklp(:,:,1:mn)-pkl(:,:,1:mn),[numel(pklp(:,:,1:mn)) 1 1]));
            %             drelnorm=dnormn/norm(reshape(pklp(:,:,1:mn),[numel(pklp(:,:,1:mn)) 1 1]));
            %complete norm
            dnormn=norm(pj(1:Nom1*Nom2*n+np,1)-pjp(1:Nom1*Nom2*n+np,1));%+norm(fi);
            drelnorm=dnormn/norm(pjp(1:Nom1*Nom2*n+np,1));
%             [dnormn drelnorm norm(fi)]
            if dnormn>dnorm
                kNRI=kNRI+1;
            end
            %display results
            dnorm=dnormn;
            %par.num.conv
            if par.num.CONV==1
                if np==1
                    conver(end+1,:)=[lcont k_nit kNRI NBD NOS par.system.(par.parlist{cpari(1)}) nan pjp(Nom1*Nom2*n+1) dnorm drelnorm];
                elseif np==2
                    conver(end+1,:)=[lcont k_nit kNRI NBD NOS par.system.(par.parlist{cpari(1)}) par.system.(par.parlist{cpari(2)}) pjp(Nom1*Nom2*n+1) dnorm drelnorm];
                else
                    conver(end+1,:)=[lcont k_nit kNRI NBD NOS  pjp(Nom1*Nom2*n+1) pjp(Nom1*Nom2*n+2) dnorm drelnorm];
                end
            end
            %calculate time for new step
            if newstep && ~isempty(pcsten)
                Ditst=pcsten;
            end
            if par.num.disp(1)==1 && par.num.disp(2)==1
                if np==1
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9 num2str(par.system.(par.parlist{cpari(1)}),'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(dnorm,'%10.5e\n') 9 9 num2str(drelnorm,'%10.5e\n') 9 9 num2str(norm(fi),'%10.5e\n')]);
                elseif np==2
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9 num2str(par.system.(par.parlist{cpari(1)}),'%10.5e\n') 9 num2str(par.system.(par.parlist{cpari(2)}),'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(dnorm,'%10.5e\n') 9 9 num2str(drelnorm,'%10.5e\n') 9 9 num2str(norm(fi),'%10.5e\n')]);
                else
                    disp([int2str(lcont) 9 num2str(Ditst,5) 9 int2str(k_nit) 9 int2str(kNRI) 9 9 int2str(NBD) 9 9 int2str(NOS) 9 9  9  9 num2str(pjp(Nom1*Nom2*n+1)/2/pi,'%10.5e\n') 9 num2str(pjp(Nom1*Nom2*n+2)/2/pi,'%10.5e\n') 9 9 num2str(dnorm,'%10.5e\n') 9 9 num2str(drelnorm,'%10.5e\n') 9 9 num2str(norm(fi),'%10.5e\n')]);
                end
                %                 minths
                PL=PL+1;
            end
            if newstep && ~isempty(pcsten);pcsten=[];newstep=false;end
            %condition
            if PL>NPL
                PL=0;
            end
            %end of iteration
            pitsten=toc(pitstst);
            %next step
            k_nit=k_nit+1;
            cpoint=true;
        end
    else
        %rebuild first Jacobian for continuation
        J=spalloc(Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+np,ceil(0.1*(Nom1*Nom2*n+1)^2));
        %define sparse array for Jacobians
        if ~knpisdelayk
            J(1:Nom1*Nom2*n+2,1:Nom1*Nom2*n+2)=jac_collocation_2dtori(u0j,par);
            %include continuation conditions (only if this is not the first correction)
            
            %derivatives w.r.t. parameters
            J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+1:Nom1*Nom2*n+2+np)=jac_cpari_ex(u0j,par,cpari,[],numpars);
        else
            [J(1:Nom1*Nom2*n+2,1:Nom1*Nom2*n+2),outtau,outp]=jac_collocation_2dtori(u0j,par);
            %include continuation conditions (only if this is not the first correction)
            out{1}=outtau;
            out{2}=outp;
            
            %derivatives w.r.t. parameters
            J(1:Nom1*Nom2*n+2+np,Nom1*Nom2*n+2+1:Nom1*Nom2*n+2+np)=jac_cpari_ex(u0j,par,cpari,out,numpars);
        end
        %set condition (%probably error here for two param continuation!
        if np==2
            J(Nom1*Nom2*n+2+1,1:Nom1*Nom2*n+2+np)=par.system.sys_cond(u0j,par,cpari,numpars,out);
        end
        %set distance condition for pseudo arclength (upj: previous solution in the branch)
        J(Nom1*Nom2*n+2+np,1:Nom1*Nom2*n+2+np)=2*(u0j'-upj');
        dnorm=0;
        brinit=false;
        cpoint=false;
    end
    if exist('cpari','var')
        np=length(cpari);
    else
        np=0;
    end
    %NR solution
    uj=pjp;
    %% find the direction in parameter space for continuation (prediction)
    % test the quality of the convergency
    if np>0
        if brinit || (dnorm<par.num.absTOL && par.system.(par.parlist{cpari(1)})<=par.bounds(cpari(1),2) && par.system.(par.parlist{cpari(1)})>=par.bounds(cpari(1),1) && sum(isnan(pjp))==0 && sum(isinf(pjp))==0)
            if isempty(iwi)
                %check healthy continuation
                if steptype=='co'
                    if NHC==par.num.limNH
                        NBD=0;
                        if NOS>0
                            par.num.ds=dsBD;
                        end
                        NOS=0;
                    end
                    NHC=NHC+1;
                end
                %changed allowed NRI
                aNRI=par.num.maxNRIc;
                %changed allowed NR
                maxNR=par.num.maxNRc;
            end
            % build the Jacobian using previous solutions
            if lcont==0
                %set branch property
                branch.contpar=cpari;
                %set mesh
                branch.mesh.m1=mth1;
                branch.mesh.m2=mth2;
                if isempty(iwi)
                    %extend or reset vector for continuation
                    for kp=1:np
                        uj(Nom1*Nom2*n+2+kp)=par.system.(par.parlist{cpari(kp)});
                    end
                    numpars.np=np;
                    resJac=jac_cpari_ex(uj,par,cpari,[],numpars);
                    J(:,Nom1*Nom2*n+2+1:Nom1*Nom2*n+2+np)=resJac(1:Nom1*Nom2*n+2+np-1,:);
                end
            else
                J=J(1:end-1,:);%spalloc(Nom1*Nom2*n+1,Nom1*Nom2*n+1+np,ceil(0.1*(Nom1*Nom2*n+1)^2));
            end
            
            %store the converged result
            plb(1,lcont+1)=max(abs(uj(1:n:Nom1*Nom2*n)));
            plb(2:2+np,lcont+1)=uj(Nom1*Nom2*n+2:Nom1*Nom2*n+2+np);
            % store the first converged solution
            if isempty(iwi)
                if cpoint
                    branch.points(lcont+1).profile=j2kl(uj(1:Nom1*Nom2*n));
                    branch.points(lcont+1).omega1=uj(Nom1*Nom2*n+1);
                    branch.points(lcont+1).omega2=uj(Nom1*Nom2*n+2);
                    branch.points(lcont+1).par=par;
                    if par.num.CONV==1
                        branch.points(lcont+1).conv=conver;
                    end
                end       
                %fixomega mask
                if par.num.fixom==0
                    itrunc=1:Nom1*Nom2*n+2+np;
                else
                    %remove one of the phase condition
                    itrunc=[1:Nom1*Nom2*n+par.num.fixom-1 Nom1*Nom2*n+par.num.fixom+1:Nom1*Nom2*n+2+np];
                end
                wi=zeros(Nom1*Nom2*n+2+np,1);
                wi(1)=1;
                %solve the equation
                wi(itrunc(2:end),1)=-J(itrunc(1:end-1),itrunc(2:end))\J(itrunc(1:end-1),itrunc(1));
                %normalization
                wi=wi/norm(wi);
                %initial step for excelling bifurcation
                instepEB=false;
            else
                wi=iwi;
                iwi=[];
                instepEB=true;
            end
            %if this is the first point
            if lcont==0
                for kp=1:np
                    uj(Nom1*Nom2*n+2+kp)=par.system.(par.parlist{cpari(kp)});
                end
            elseif par.num.PR && ~isempty(wpi)
                %check projection to the old tangent
                if wpi.'*wi<0
                    wi=-wi;
                end
            end
            %new prediction, TODO: maybe projection is needed
            boundaryhit=true;
            %counter to prevent infinite loop
            inflc=0;
            while boundaryhit
                if ~instepEB
                    ujp=uj+wi*par.num.ds;%/norm(wi(end-np+1:end));
                else
%                     ujp=uj+wi*par.num.Ids;
                    ujp(1:Nom1*Nom2*n+2,1)=uj(1:Nom1*Nom2*n+2,1)+wi(1:Nom1*Nom2*n+2,1)*par.num.Ids;
                    ujp(Nom1*Nom2*n+2+1,1)=uj(Nom1*Nom2*n+2+1,1)+wi(Nom1*Nom2*n+2+1,1);
                    uj(Nom1*Nom2*n+2+1,1)=ujp(Nom1*Nom2*n+2+1,1);
                end
                boundaryhit=false;
                for kp=1:np
                    if sum(sum(isinf(wi),1),2)>0 || sum(sum(isnan(wi),1),2)>0 || (ujp(Nom1*Nom2*n+2+kp)>=par.bounds(cpari(kp),1) && ujp(Nom1*Nom2*n+2+kp)<=par.bounds(cpari(kp),2))
                        par.system.(par.parlist{cpari(kp)})=ujp(Nom1*Nom2*n+2+kp);
                    else
                        if inflc<10
                            boundaryhit=true;
                        else
                            boundaryhit=false;
                        end
                        inflc=inflc+1;
                        par.num.ds=par.num.ds/2;
                        ['boundary hit on parameter ' par.parlist{cpari(kp)}]
                    end
                end
                
            end
            %set previous convergent solution
            upj=uj;
            u0j=ujp;
%             size(u0j)
            if par.num.plot(1) && ~isempty(wpi)
                set(0,'CurrentFigure',fig2);
                hold on;
                %prediction
                if length(cpari)<=1 && PNP==1
                    plot([plb(3,lcont+1) ujp(Nom1*Nom2*n+2+1)], [plb(1,lcont+1) max(abs(ujp(1:n:Nom1*Nom2*n)))],prec)
                    %plot the last continued branch part
                    if lcont>0 && cpoint
                        plot(plb(3,lcont:lcont+1),plb(1,lcont:lcont+1),'k.-');
                    end
                    %plot converged circle
                    plot(plb(3,lcont+1),plb(1,lcont+1),'go');
                else
                    %f
                    plot([plb(4,lcont+1) ujp(Nom1*Nom2*n+2+2)], [plb(3,lcont+1) ujp(Nom1*Nom2*n+1+2)],prec)
                    %plot the last continued branch part
                    if lcont>0
                        plot(plb(4,lcont:lcont+1),plb(3,lcont:lcont+1),'k.-');
                    end
                    %plot converged circle
                    plot(plb(4,lcont+1),plb(3,lcont+1),'go');
                end
                drawnow;
            end
            %store old tangent
            wpi=wi;
            %change step size
            if k_nit-1<=2 && k_nit-1>=1
                %grow stepsize
                par.num.ds=par.num.gds*par.num.ds;
            elseif k_nit-1>=3%par.num.maxNR-1
                %decrese step size
                par.num.ds=par.num.ds/par.num.gds;
            end
            %store valid step length for break down
            if (NOS==0 && NBD==0 && NHC<par.num.limNH) || NHC>=par.num.limNH
                dsBD=par.num.ds;
            end
            %increae step
            lcont=lcont+1;
            %set PL==0, change type gicing heading
            if lcont==1
                PL=0;
            end
            %steptype
            steptype='co';
            %save backup
            bkpbranch=branch;
            if ~isempty(stepsavename)
                save([stepsavename '.mat'],'bkpbranch');
            end
            brinit=false;
        else
            if lcont>0
                %store the old point
                ujold=ujp;
                %wrong point, not convergent, use the half distance, and set ds half
                if NBD<par.num.maxBRD
                    %set half ds
                    par.num.ds=par.num.ds/2;
                    %count break downs
                    NBD=NBD+1;
                    %steptype
                    steptype='db';
                else
                    %overshoot
                    if NBD==par.num.maxBRD
                        par.num.ds=dsBD*(2+NOS)*3/pi;
                    else
                        par.num.ds=2*par.num.ds;
                    end
                    NOS=NOS+1;
                    %steptype
                    steptype='do';
                end
                %new predicted point
                ujp=upj+wi*par.num.ds;%/norm(wi(end-np+1:end));
                for kp=1:np
                    par.system.(par.parlist{cpari(kp)})=ujp(Nom1*Nom2*n+1+kp);
                end
                u0j=ujp;
                if par.num.plot(1)
                    set(0,'CurrentFigure',fig2);
                    hold on;
                    %prediction
                    if np<=1 && PNP==1
                        plot([plb(3,lcont) ujp(Nom1*Nom2*n+2+1)], [plb(1,lcont) max(abs(ujp(1:n:Nom1*Nom2*n)))],prec)
                        %plot the unconverged cicle
                        plot(ujold(Nom1*Nom2*n+2+1),max(abs(ujold(1:n:Nom1*Nom2*n))),'ro');
                    else
                        plot([plb(4,lcont) ujp(Nom1*Nom2*n+1+3)],[plb(3,lcont) ujp(Nom1*Nom2*n+2+1)],prec);
                        %plot the unconverged cicle
                        plot(ujold(Nom1*Nom2*n+1+3),ujold(Nom1*Nom2*n+2+1),'ro');
                    end
                    drawnow;
                end
                
            else
                if nargout==1
                    error('The initial condition is not convergent!') ;
                else
                    success=false;
                end
            end
        end
    elseif ~(((k_nit<=par.num.maxNR && kNRI<=par.num.maxNRI && lcont==0) || (k_nit<=par.num.maxNRc && kNRI<=par.num.maxNRIc && lcont>0)) && dnorm<par.num.absTOL)
        if nargout==1
            error('The initial condition is not convergent!') ;
        else
            success=false;
            %store the solution for a point
            branch.points(lcont+1).profile=j2kl(uj(1:Nom1*Nom2*n));
            branch.points(lcont+1).omega1=uj(Nom1*Nom2*n+1);
            branch.points(lcont+1).omega2=uj(Nom1*Nom2*n+2);
            branch.points(lcont+1).par=par;
            if par.num.CONV==1
                branch.points(lcont+1).conv=conver;
            end
            %a step to end the outest while
            lcont=lcont+1;
        end
    else
        %store the solution for a point
        branch.points(lcont+1).profile=j2kl(uj(1:Nom1*Nom2*n));
        branch.points(lcont+1).omega1=uj(Nom1*Nom2*n+1);
        branch.points(lcont+1).omega2=uj(Nom1*Nom2*n+2);
        branch.points(lcont+1).par=par;
        if par.num.CONV==1
            branch.points(lcont+1).conv=conver;
        end
        %a step to end the outest while
        lcont=lcont+1;

        success=true;
    end
            newstep=true;
    pcsten=toc(pcstst);
end
toc;
if NOS>=par.num.maxOVS
    disp('Continuation has stopped after break downs and overshoot trials!');
end
beep;
return
%% NESTED FUNCTIONS

    function uj=kl2j(varargin)
        if length(varargin)==1
            n1=size(varargin{1},1);
            n2=size(varargin{1},2);
        elseif length(varargin)==2
            n1=varargin{2};
            n2=n1;
        else
            n1=varargin{2};
            n2=varargin{3};
        end
        pukl=varargin{1};
        uj=reshape(permute(pukl,[3 1 2]),[n1*n2*size(pukl,3) 1]);
    end

    function pukl=j2kl(varargin)
        if length(varargin)==1
            pukl=permute(reshape(varargin{1},[n, Nom1, Nom2]),[2 3 1]);
        elseif length(varargin)==2
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{2}]),[2 3 1]);
        elseif length(varargin)==3
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{3}]),[2 3 1]);
        elseif length(varargin)==4
            pukl=permute(reshape(varargin{1},[varargin{4}, varargin{2}, varargin{3}]),[2 3 1]);
        end
    end

   
%Lagrange polynomial
    function res=lag(eps,p)
        ms=0:P;
        ms=ms(ms~=p);
        res=zeros(1,length(eps));
        for klag=1:length(eps)
            res(1,klag)=prod((P*eps(klag)-ms)./(p-ms));
        end
    end

%derivative of Lagrange polynomial
    function res=dlag(eps,p)
        res=zeros(1,length(eps));
        ns=0:P;
        ns=ns(ns~=p);
        %         res=0;
        for klag=1:length(eps)
            for kd=1:P
                ms=ns(ns~=ns(kd));
                res(1,klag)=res(1,klag)+P/(p-ns(kd))*prod((P*eps(klag)-ms)./(p-ms));
            end
        end
    end


    function Jres=jac_cond(vj,cpari)
        %REMARK: analyticalm
        %TODO: it is not finished for multi delays
        %determine minimum
        [~,klmins,chminkl,thmin12,himin,gimin,phiimin]=fun_condition(vj,par);
        %number of teeth
        Z=par.system.Z;
        %collocation interval length
        dth1=2*pi/NC1;
        %collocation interval
        dth2=2*pi/NC2;
        %select the minimum good chip thickness
        [~,mI]=min(himin+(2-gimin));
        %         mI=1;
        %declare lagrange at tmin1 tmin2
        ki=mI;
        ki=2;
        %         for ki=1:Z
        %product of other minima
        promin=1;%prod(himin([1:ki-1 ki+1:Z]))*1000^(Z-1);
        %normal at the minimum time instant
        nimin=[sin(phiimin(ki))*sin(par.system.kappa);sin(par.system.kappa)*cos(phiimin(ki));-cos(par.system.kappa);];
        %coefficient of the derivative
        cderi=permute([par.system.us zeros(size(par.system.us))].'*nimin,[2 3 1])/par.system.fZ;
        tm1lag=zeros(P+1,1);
        tm2lag=zeros(P+1,1);
        %mod minimu
        thmin12_mod=[mod(thmin12(ki,1),dth1)/dth1 mod(thmin12(ki,2),dth2)/dth2];
        %lagrange approximation at the minimum
        for KP=0:P
            tm1lag(KP+1)=lag(thmin12_mod(1),KP);
            tm2lag(KP+1)=lag(thmin12_mod(2),KP);
        end
        %shape functions at tmin1 and tmin2
        tm12Ppq=tm1lag*tm2lag.';
        %preallocate
        Jres=spalloc(1,Nom1*Nom2*n+2+np,ceil(0.1*(Nom1*Nom2*n+1)^2));
        Jres(1,1:Nom1*Nom2*n+1)=Jres(1,1:Nom1*Nom2*n+1)+sparse(ones((P+1)^2*n+1,1),[IR(ceil((mod(klmins(ki,1)-1,NC1*P)+1)/P),ceil((mod(klmins(ki,2)-1,NC2*P)+1)/P));Nom1*Nom2*n+1],[reshape(promin*repmat(cderi,[P+1 P+1 1]).*repmat(tm12Ppq,[1 1 n]),[(P+1)^2*n 1]);0;]);
        %frequencies
        om1=vj(Nom1*Nom2*n+1);
        om2=vj(Nom1*Nom2*n+2);
        %calculate the delay, condit
        tau_pc=par.system.sys_tau(vj,par);
        %delay tmin
        thmin12_tau=[mod(thmin12(ki,1)-om1*tau_pc,2*pi) mod(thmin12(ki,2)-om2*tau_pc,2*pi)];
        %mod minimu
        thmin12_tau_mod=[mod(thmin12_tau(1),dth1)/dth1 mod(thmin12_tau(2),dth2)/dth2];
        %declare lagrange at tmin1 tmin2, tau
        tm1lag_tau=zeros(P+1,1);
        tm2lag_tau=zeros(P+1,1);
        dtm2lag_tau_dtm2=zeros(P+1,1);
        for KP=0:P
            tm1lag_tau(KP+1)=lag(thmin12_tau_mod(1),KP);
            tm2lag_tau(KP+1)=lag(thmin12_tau_mod(2),KP);
            dtm2lag_tau_dtm2(KP+1)=dlag(thmin12_tau_mod(2),KP);
        end
        %shape functions at tmin1 and tmin2
        tm12Ppq_tau=tm1lag_tau*tm2lag_tau.';
        %shape functions (dth2) at tmin1 and tmin2
        dtm12Ppq_tau_dtm2=tm1lag_tau*dtm2lag_tau_dtm2.';
        %delayed interval indeces
        taukl=ceil(thmin12_tau/dth1);
        %store in the Jacobian for state variables
        Jres(1,1:Nom1*Nom2*n+1)=Jres(1,1:Nom1*Nom2*n+1)+sparse(ones((P+1)^2*n+1,1),[IR(taukl(1),taukl(2));Nom1*Nom2*n+1;],[reshape(-promin*repmat(cderi,[P+1 P+1 1]).*repmat(tm12Ppq_tau,[1 1 n]),[(P+1)^2*n 1]);0;]);
        %state
        vkl=j2kl(vj(1:Nom1*Nom2*n));
        %delayed segment
        vkltau=vkl((taukl(1)-1)*P+(1:P+1),(taukl(2)-1)*P+(1:P+1),:);
        %qtauimin
        qtauimimn=sum(sum(vkltau.*repmat(dtm12Ppq_tau_dtm2,[1 1 n]),1),2);
        %frequency dependency, \omega_2
        Jres(1,Nom1*Nom2*n+2)=Jres(1,Nom1*Nom2*n+2)-permute(cderi,[1 3 2])*permute(qtauimimn,[3 1 2])*(-tau_pc);
        %         end
    end
%linear indexes for representation vector at (k,l) interval
    function res=IR(k,l)
        res=reshape(repmat(((l-1)*P+(1:P+1)-1)*(NC1*P+1)*n,[P+1 1 n])+repmat(permute(((k-1)*P+(1:P+1)-1)*n,[2 1 3]),[1 P+1 n])+repmat(permute(1:n,[1 3 2]),[P+1 P+1 1]),[(P+1)^2*n 1]);
    end

%linear indexes for representation vector at (k,l) interval shifted by r1 and r2
    function res=IRs(k,l,r1,r2)
        res=reshape(repmat((mod((l-1)*P+(1:P+1)-r2-1,NC1*P+1)+1-1)*(NC1*P+1)*n,[P+1 1 n])+repmat(permute((mod((k-1)*P+(1:P+1)-1-r1,NC2*P+1)+1-1)*n,[2 1 3]),[1 P+1 n])+repmat(permute(1:n,[1 3 2]),[P+1 P+1 1]),[(P+1)^2*n 1]);
    end

end

