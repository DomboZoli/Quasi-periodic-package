function [Jc,varargout]=jac_collocation_2dtori(pj,par,varargin)
% build collocational Jacobian for quasiperiodic vibration
%% VERSION v10.1
%v2.0: not brute force building of collocation Jacobianr
%v3.0: working segmented non brute code
%v4.0: paralell computing for the segmented code%
%v4.1: paralell computing for the segmented code ineer loop for smaller storage
%v4.2: two unknown freed frequencies
%v4.2b: v4.2 modified for 1dof turning
%v5.0b: analytical derivatives, where it is possible
%v5.1b: one polynom should be enough
%v6.0: one polynom should be enough
%v6.1: all derivatives are replaced
%v6.2: parfor kJ
%v6.3: smart sparse
%v7.0: for normal stability it serves entire jacobian, separated from the delays
%v7.1: multiple constant delays
%v7.2: milling with NFOs
%v7.3: sysderi modified
%v7.4: only real vibration frequency
%v8.0: introduction of sys_rhs and sys_deri
%v8.1: selectable scheme, par.num.scheme, 1: Chebishev, 2: Legendre, 3:Lobatto
%v8.2: corrected frequency derivatives
%v8.3: corrected bug in frequency derivatives
%v9.0: idealistic case distinguished from flyover, par.system.NFO==0: ideal, par.system.NFO>0 flyover
%v9.1: dim and dof in the par.system list, def_taus introduced
%v9.2: save G for cpari
%v9.3: debugged
%v9.4: phase condition over kl instead of ij
%v9.5: phase condition over cij recap
%v10.0: finalized parallellization
%v10.1: debug 1
%% INPUT
%   pj: actual (jth) iteration step, [puij; pom;]
%       puij: (p: not iterated) representation points
%       pom: (p: not iterated) vibration frequency
%   par: parameterlist (see parameterlist file)
%   varargin{1}=<stab>: jacobian for stability without the boundary conditions, and separated present and delayed state if stab==1
%% OUTPUT
%   Jc: colocation Jacobian
%   varargout{1}=Jctau: delayed colocation Jacobian
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% INCONSTANTS
%set up phase condition reference solution
uj0=par.num.uj0;
%position tolerance
uTOL=1e-8;
%brute force to build up the jacobian
brute=false;
%check if Jc for NR or normal stability
if ~isempty(varargin)
    stab=varargin{1};
else
    %for NR
    stab=0;
end
%finite difference constants
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
%collocation order
P=par.num.p;
cpps=zeros(7,7);
wpps=zeros(7,7);
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
        % Chebyshev nodes of the first kind
        for kch=1:7;cpps(kch,1:kch)=flip((1+cos((2*(1:kch)-1)/(2*kch)*pi))/2);end
        % Chebyshev weights for the gauss quadrature, sqrt and /2 is due to coordinate transform from (-1,1) to (0,1)
        for kch=1:7;wpps(kch,1:kch)=ones(1,kch).*sqrt(1-(2*(cpps(kch,1:kch)-1/2)).^2).*pi/kch/2;end
    case 2
        %Legendre
        cpps(1,1)=0.5;
        cpps(2,1:2)=[0.211325, 0.788675];
        cpps(3,1:3)=[0.112702, 0.5, 0.887298];
        cpps(4,1:4)=[0.0694318, 0.330009, 0.669991, 0.930568];
        cpps(5,1:5)=[0.0469101, 0.230765, 0.5, 0.769235, 0.95309];
        cpps(6,1:6)=[0.0337652, 0.169395, 0.38069, 0.61931, 0.830605, 0.966235];
        cpps(7,1:7)=[0.025446, 0.129234, 0.297077, 0.5, 0.702923, 0.870766, 0.974554];
        %weights
        wpps(1,1)=2;
        wpps(2,1:2)=[1 1];
        wpps(3,1:3)=[5/9 8/9 5/9];
        wpps(4,1:4)=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
        wpps(5,1:5)=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
        %due to coordinate transform from (-1,1) to (0,1)
        wpps=wpps/2;
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
%quadrature
wp=wpps(P,1:P);
%Legendre
% switch P

%% parameters
%finite difference constants
m=par.num.dp;
%order of the cental difference
mp=round(m/2);
%collocation intervals for 2D tori both directions
NC1=par.num.NC1;
NC2=par.num.NC2;
%number of representation points
Nom1=NC1*P+1;
Nom2=NC2*P+1;
%collocation interval length
dth1=2*pi/NC1;
%collocation interval
dth2=2*pi/NC2;
%system DOF
dof=par.system.dof;
%dimension
n=2*dof;
%representation stepss
h1=2*pi/(NC1*P);
h2=2*pi/(NC2*P);
%all shape functions in all collocation points
cPpq=zeros(P,P,P+1,P+1);
cd1Ppq=zeros(P,P,P+1,P+1);
cd2Ppq=zeros(P,P,P+1,P+1);
td1Ppq=zeros(P+1,P+1,P+1,P+1);
td2Ppq=zeros(P+1,P+1,P+1,P+1);
%interval points
th1i=((1:NC1+1)-1)*dth1;
th2j=((1:NC2+1)-1)*dth2;
%collocation points
c1ip=reshape(repmat(th1i(1:end-1),[P,1])+repmat(cp'*dth1,[1 NC1]),[NC1*P,1]).';
c2jp=reshape(repmat(th2j(1:end-1),[P,1])+repmat(cp'*dth2,[1 NC2]),[NC2*P,1]).';
%check weather tau is in continuation parameters
%calculate the delays
[~,npt]=par.system.sys_tau([],par);
%delay parameter index
dpi=find(npt==par.num.cpari,1);
%the parameter is a delay? if yes G must be stored for cpari
knpisdelay=~isempty(dpi) || numel(par.num.cpari)>1;

%% main code
%preallocation
if stab
    Jc=spalloc(NC1*NC2*P^2*n,Nom1*Nom2*n,2*Nom1*Nom2*P*n);
    Jctau=spalloc(NC1*NC2*P^2*n,((Nom1+1)*(Nom2+1))*n,2*Nom1*Nom2*P*n);
else
    Jc=spalloc(Nom1*Nom2*n+2,Nom1*Nom2*n+2,2*Nom1*Nom2*P*n);
end
%phase condition shape function matrix
if par.num.pctype~=0
    Mpresder1=spalloc((NC1*P+1)*(NC2*P+1)*n+1,NC1*P*NC2*P*n+1,Nom1*Nom2*P*n);
    Mpresder2=spalloc((NC1*P+1)*(NC2*P+1)*n+1,NC1*P*NC2*P*n+1,Nom1*Nom2*P*n);
end
%brute force derivation by central differenc schemes
% for kk=1:Nom1*Nom2*n
if brute
    parfor kk=1:Nom1*Nom2*n+2
        %     for kk=1:Nom1*Nom2*n
        %one coloumn
        Jccol=zeros(Nom1*Nom2*n+2,1);
        for km=-mp:mp
            Jccol=Jccol+etas(km+mp+1)*fun_invariance_2dtori(pj(1:Nom1*Nom2*n+2)+[zeros((kk-1),1); km; zeros(Nom1*Nom2*n+2-kk,1)]*uTOL,par);
        end
        %     [kk Nom1*Nom2*n]
        Jc(:,kk)=Jccol/uTOL;
    end
    %calculate the delays,
    taus=par.system.sys_tau([],par);
    %number of delays
    Ntau=numel(taus);
    %paralell computing variables
    Ppc=P;
    NC1pc=NC1;
    NC2pc=NC2;
    %invariable
    npc=n;
    %produce the polinomial series at collocational points
    clag=zeros(P+1,P);
    cdlag=zeros(P+1,P);
    
    for kP=0:P
        clag(kP+1,:)=lag(cp,kP);
        cdlag(kP+1,:)=dlag(cp,kP);
    end
    for pp=1:P+1
        for qq=1:P+1
            cPpq(:,:,pp,qq)=clag(pp,:)'*clag(qq,:);
            cd1Ppq(:,:,pp,qq)=cdlag(pp,:)'*clag(qq,:)/dth1;
            cd2Ppq(:,:,pp,qq)=clag(pp,:)'*cdlag(qq,:)/dth2;
        end
    end
    cPpqpc=cPpq;
    %data for external functions
    numpars.Nom1=Nom1;
    numpars.Nom2=Nom2;
    numpars.P=P;
    numpars.n=n;
    numpars.cPpq=cPpq;
    numpars.cd1Ppq=cd1Ppq;
    numpars.cd2Ppq=cd2Ppq;
    %vibration frequency
    om1=real(pj(Nom1*Nom2*n+1,1));
    om2=real(pj(Nom1*Nom2*n+2,1));
else
    %produce the polinomial series at collocational points
    clag=zeros(P+1,P);
    cdlag=zeros(P+1,P);
    %produce the polinomial series at representation points
    tlag=zeros(P+1,P+1);
    tdlag=zeros(P+1,P+1);
    tp=((1:P+1)-1)/P;
    
    for kP=0:P
        clag(kP+1,:)=lag(cp,kP);
        cdlag(kP+1,:)=dlag(cp,kP);
        tlag(kP+1,:)=lag(tp,kP);
        tdlag(kP+1,:)=dlag(tp,kP);
    end
    % collocation interval length
    dth1=2*pi/NC1;
    dth2=2*pi/NC2;
    for pp=1:P+1
        for qq=1:P+1
            cPpq(:,:,pp,qq)=clag(pp,:)'*clag(qq,:);
            cd1Ppq(:,:,pp,qq)=cdlag(pp,:)'*clag(qq,:)/dth1;
            cd2Ppq(:,:,pp,qq)=clag(pp,:)'*cdlag(qq,:)/dth2;
            td1Ppq(:,:,pp,qq)=tdlag(pp,:)'*tlag(qq,:)/dth1;
            td2Ppq(:,:,pp,qq)=tlag(pp,:)'*tdlag(qq,:)/dth2;
        end
    end
    cPpqpc=cPpq;
    %calculate the delay, TODO: for general delays
    taus=par.system.sys_tau([],par);
    %number of tau
    Ntau=numel(taus);
    %vibration frequency
    om1=real(pj(Nom1*Nom2*n+1,1));
    om2=real(pj(Nom1*Nom2*n+2,1));
    %derive pkl
    pkl=j2kl(pj(1:Nom1*Nom2*n));
    %system matrices
    n=2*dof;
    %paralell computing variables
    Ppc=P;
    NC1pc=NC1;
    NC2pc=NC2;
    %data for external functions
    numpars.Nom1=Nom1;
    numpars.Nom2=Nom2;
    numpars.P=P;
    numpars.n=n;
    numpars.cPpq=cPpq;
    numpars.cd1Ppq=cd1Ppq;
    numpars.cd2Ppq=cd2Ppq;
    %invariable
    npc=n;
    %\om_t \pd/\pd\theta1
    %present states 6D
    prestate=zeros(P^2*n^2*(P+1)^2,NC1,NC2);
    %present derivative matrix
    if par.num.pctype~=0
        presder1=zeros(P^2*n^2*(P+1)^2,NC1,NC2);
        presder2=zeros(P^2*n^2*(P+1)^2,NC1,NC2);
    end
    %delayed state
    %     delstate=zeros(P,P,n,P+2,P+2,n,NC1*NC2);
    delstate=zeros(Ppc^2*npc^2*(Ppc+1)^2,NC1,NC2,Ntau);
    %store G if tau is continuation parameter
    if knpisdelay
        SGijkl=zeros(Ppc,Ppc,npc,NC1,NC2,npc,1);
        SGijkltau=zeros(Ppc,Ppc,npc,NC1,NC2,npc,Ntau);
    end
    %declare indeces
    i_1=zeros(1,NC1,Ntau);
    i_2=zeros(1,NC2,Ntau);
    Dc1=zeros(1,NC1,Ntau);
    Dc2=zeros(1,NC2,Ntau);
    for ktau=1:Ntau
        for kJ=1:NC1
            %delayed mesh time (checked!)
            cijtau_1=mod((kJ-1)*dth1-om1*taus(ktau),2*pi);
            %calculate minimum modulo to a coll point
            i_1(:,kJ,ktau)=floor(cijtau_1/h1);
            %difference bitween the
            Dc1(:,kJ,ktau)=cijtau_1/h1-i_1(:,kJ,ktau);
            %if it is about stability calculation
            if stab
                i_1(:,kJ,ktau)=(kJ-1)*P;
                %shift in second coordinate
                stau(1,ktau)=ceil(om1*taus(ktau)/h1);
            end
        end
        for lJ=1:NC2
            %delayed mesh time (checked!)
            cijtau_2=mod((lJ-1)*dth2-om2*taus(ktau),2*pi);
            %calculate minimum modulo to a coll point
            i_2(:,lJ,ktau)=floor(cijtau_2/h2);
            %difference bitween the
            Dc2(:,lJ,ktau)=cijtau_2/h2-i_2(:,lJ,ktau);
            %if it is about stability calculation
            if stab
                i_2(:,lJ)=(lJ-1)*P;
                %shift in second coordinate
                stau(2,ktau)=ceil(om2*taus(ktau)/h2);
            end
        end
    end
    %computing segments
    %PARFOR
    %                 for kJ=1:NC1
    parfor kJ=1:NC1
        %principle period frequency (known, till there is no flyover this is the tooth passing frequency)
        
        for lJ=1:NC2
            %cut out (k,l) (Ppc+1)^2*npc interval
            pklJ=pkl((kJ-1)*Ppc+1:kJ*Ppc+1,(lJ-1)*Ppc+1:lJ*Ppc+1,:);
            %system at (kJ,lJ) interval
            omdth=zeros(Ppc,Ppc,npc,Ppc+1,Ppc+1,npc);
            for kn=1:npc
                omdth(:,:,kn,:,:,kn)=om1*cd1Ppq+om2*cd2Ppq;
            end
            %declare force
            Gijkl=zeros(Ppc,Ppc,npc,1,1,npc);
            Gijkltau=zeros(Ppc,Ppc,npc,1,1,npc,Ntau);
            %present state
            pcijJ=cpstate_ex(pklJ,numpars,0);
            pcijJtauom=zeros(Ppc,Ppc,npc,Ntau);
            %delayed state
            for ktau=1:Ntau
                pcijJtauom(:,:,:,ktau)=cptaustate_ex(pkl,om1*taus(ktau),om2*taus(ktau),(kJ-1)*Ppc+(1:(Ppc+1)),(lJ-1)*Ppc+(1:(Ppc+1)),numpars,0);
            end
            for kF=1:Ppc
                for lF=1:Ppc
                    Gijkl(kF,lF,:,1,1,:)=permute(par.system.sys_deri(c1ip((kJ-1)*Ppc+kF),c2jp((lJ-1)*Ppc+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1),[3,4,1,5,6,2]);
                    for ktau=1:Ntau
                        Gijkltau(kF,lF,:,1,1,:,ktau)=permute(par.system.sys_deri(c1ip((kJ-1)*Ppc+kF),c2jp((lJ-1)*Ppc+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1+ktau),[3,4,1,5,6,2]);
                    end
                end
            end
            if knpisdelay
                SGijkl(:,:,:,kJ,lJ,:,1)=Gijkl;
                SGijkltau(:,:,:,kJ,lJ,:,:)=Gijkltau;
            end
            %present Jacobian, inlude partial derivatices, 6D
            sysijkl=omdth-(repmat(Gijkl,[1,1,1,Ppc+1,Ppc+1,1])).*repmat(permute(cPpqpc,[1 2 5 3 4 6]),[1 1 npc 1 1 npc]);
            %store
            prestate(:,kJ,lJ)=reshape(sysijkl,[Ppc^2*npc^2*(Ppc+1)^2,1]);
            if par.num.pctype~=0
                presder1(:,kJ,lJ)=reshape(repmat(permute(cd1Ppq,[1 2 5 3 4 6]),[1 1 npc 1 1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]);
                presder2(:,kJ,lJ)=reshape(repmat(permute(cd2Ppq,[1 2 5 3 4 6]),[1 1 npc 1 1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]);
            end
            sysijkl_tau=zeros(Ppc,Ppc,npc,Ppc+1,Ppc+1,npc,Ntau);
            for ktau=1:Ntau
                sysijkl_tau(:,:,:,:,:,:,ktau)=-repmat(Gijkltau(:,:,:,:,:,:,ktau),[1,1,1,Ppc+1,Ppc+1,1]).*repmat(permute(cPpqpc,[1,2,5,3,4,6]),[1,1,npc,1,1,npc]);
            end
            %delayed Jacobian, including derivatives, 7D
            delstate(:,kJ,lJ,:)=reshape(permute(sysijkl_tau,[1 2 3 4 5 6 7]),[Ppc^2*npc^2*(Ppc+1)^2,1,1,Ntau]);
        end
    end
    
    
    if knpisdelay
        varargout{1}=SGijkltau;
        varargout{2}=SGijkl;
    end
    %store in the Jacobian, TODO: paralellize, due to the storing parfor is not effective here
    %     parfor kJ=1:NC1pc
    
    for kJ=1:NC1pc
        if stab
            prJc=spalloc(NC1*NC2*P^2*n,Nom1*Nom2*n,2*Nom1*Nom2*P*n);
            prJctau=spalloc(NC1*NC2*P^2*n,((Nom1+1)*(Nom2+1))*n,2*Nom1*Nom2*P*n);
        else
            prJc=spalloc(Nom1*Nom2*n+2,Nom1*Nom2*n+2,2*Nom1*Nom2*P*n);
            prJctau=spalloc(Nom1*Nom2*n+2,Nom1*Nom2*n+2,2*Nom1*Nom2*P*n);
        end
        if par.num.pctype~=0
            prMpresder1=spalloc((NC1*P+1)*(NC2*P+1)*n+1,(NC1*P)*(NC2*P)*n+1,Nom1*Nom2*P*n);
            prMpresder2=spalloc((NC1*P+1)*(NC2*P+1)*n+1,(NC1*P)*(NC2*P)*n+1,Nom1*Nom2*P*n);
        end
        %due to the storing parfor is not effective here
        for lJ=1:NC2pc
            spi=zeros(P^2*n^2*(P+1)^2+1,1);
            spj=zeros(P^2*n^2*(P+1)^2+1,1);
            spji=zeros(P^2*n^2*(P+1)^2+1,1);
            spjj=zeros(P^2*n^2*(P+1)^2+1,1);
            %set last element
            spi(P^2*n^2*(P+1)^2+1)=size(prJc,1);
            spj(P^2*n^2*(P+1)^2+1)=size(prJc,2);
            %collocation indeces
            spi(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(fun_IC(kJ,lJ,numpars),[1 (P+1)^2*n]),[P^2*n^2*(P+1)^2,1]);
            %representation indeces of present state
            spj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(fun_IR(kJ,lJ,numpars).',[P^2*n 1]),[P^2*n^2*(P+1)^2,1]);
            %representation indeces of present state
            spji(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(fun_IR(kJ,lJ,numpars),[1 P^2*n]),[P^2*n^2*(P+1)^2,1]);
            spjj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(fun_IC(kJ,lJ,numpars).',[(P+1)^2*n 1]),[P^2*n^2*(P+1)^2,1]);
            %store present state in Jacobian
            prJc=prJc+sparse(spi,spj,[prestate(:,kJ,lJ); 0;]);
            %store derivative shape functions for phase conditions
            if par.num.pctype~=0
                spji(P^2*n^2*(P+1)^2+1)=size(prMpresder1,1);
                spjj(P^2*n^2*(P+1)^2+1)=size(prMpresder2,2);
                prMpresder1=prMpresder1+sparse(spji,spjj,[presder1(:,kJ,lJ); 0;]);
                prMpresder2=prMpresder2+sparse(spji,spjj,[presder2(:,kJ,lJ); 0;]);
            end
            %delayed representation indeces
            for ktau=1:Ntau
                %patch indeces
                if ~stab
                    inds1P2_w1=mod(i_1(1,kJ,ktau)+(0:P),NC1*P)+1;
                    inds1P2_w2=mod(i_1(1,kJ,ktau)+(0:P)+1,NC1*P)+1;
                    inds2P2_w1=mod(i_2(1,lJ,ktau)+(0:P),NC2*P)+1;
                    inds2P2_w2=mod(i_2(1,lJ,ktau)+(0:P)+1,NC2*P)+1;
                else
                    inds1P2_w1=i_1(1,kJ,ktau)+(0:P)+1;
                    inds1P2_w2=i_1(1,kJ,ktau)+(0:P)+1+1;
                    inds2P2_w1=i_2(1,lJ,ktau)+(0:P)+1;
                    inds2P2_w2=i_2(1,lJ,ktau)+(0:P)+1+1;
                end
                %TODO: check this part consistency, it works well but the kpC and lpC for cycles seems unnecessary or at least inefficient
                W1_1=repmat(abs(1-abs(Dc1(1,kJ,ktau))),[P P]);
                W1_2=repmat(abs(1-abs(Dc2(1,lJ,ktau))),[P P]);
                DELEYIR_w11=repmat(permute(repmat((inds2P2_w1-1)*(NC1*P+1)*n,[Ppc+1 1 n])+repmat(permute((inds1P2_w1-1)*n,[2 1 3]),[1 Ppc+1 n])+repmat(permute(1:n,[1 3 2]),[Ppc+1 Ppc+1 1]),[4 5 6 1 2 3]),[P P 1 1 1 1]);
                DELEYIR_w21=repmat(permute(repmat((inds2P2_w1-1)*(NC1*P+1)*n,[Ppc+1 1 n])+repmat(permute((inds1P2_w2-1)*n,[2 1 3]),[1 Ppc+1 n])+repmat(permute(1:n,[1 3 2]),[Ppc+1 Ppc+1 1]),[4 5 6 1 2 3]),[P P 1 1 1 1]);
                DELEYIR_w12=repmat(permute(repmat((inds2P2_w2-1)*(NC1*P+1)*n,[Ppc+1 1 n])+repmat(permute((inds1P2_w1-1)*n,[2 1 3]),[1 Ppc+1 n])+repmat(permute(1:n,[1 3 2]),[Ppc+1 Ppc+1 1]),[4 5 6 1 2 3]),[P P 1 1 1 1]);
                DELEYIR_w22=repmat(permute(repmat((inds2P2_w2-1)*(NC1*P+1)*n,[Ppc+1 1 n])+repmat(permute((inds1P2_w2-1)*n,[2 1 3]),[1 Ppc+1 n])+repmat(permute(1:n,[1 3 2]),[Ppc+1 Ppc+1 1]),[4 5 6 1 2 3]),[P P 1 1 1 1]);
                %set last element
                spi(P^2*n^2*(P+1)^2+1)=size(prJctau,1);
                spj(P^2*n^2*(P+1)^2+1)=size(prJctau,2);
                %w_11
                spj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(DELEYIR_w11,[1 1 n 1 1 1]),[P^2*n^2*(Ppc+1)^2,1]);
                prJctau=prJctau+sparse(spi,spj,[reshape(repmat(W1_1.*W1_2,[1 1 npc Ppc+1 Ppc+1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]).*delstate(:,kJ,lJ,ktau);0;]);
                %w_12
                spj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(DELEYIR_w12,[1 1 n 1 1 1]),[P^2*n^2*(Ppc+1)^2,1]);
                prJctau=prJctau+sparse(spi,spj,[reshape(repmat(W1_1.*(1-W1_2),[1 1 npc Ppc+1 Ppc+1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]).*delstate(:,kJ,lJ,ktau);0;]);
                %w_21
                spj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(DELEYIR_w21,[1 1 n 1 1 1]),[P^2*n^2*(Ppc+1)^2,1]);
                prJctau=prJctau+sparse(spi,spj,[reshape(repmat((1-W1_1).*W1_2,[1 1 npc Ppc+1 Ppc+1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]).*delstate(:,kJ,lJ,ktau);0;]);
                %w_22
                spj(1:P^2*n^2*(P+1)^2,1)=reshape(repmat(DELEYIR_w22,[1 1 n 1 1 1]),[P^2*n^2*(Ppc+1)^2,1]);
                prJctau=prJctau+sparse(spi,spj,[reshape(repmat((1-W1_1).*(1-W1_2),[1 1 npc Ppc+1 Ppc+1 npc]),[Ppc^2*npc^2*(Ppc+1)^2,1]).*delstate(:,kJ,lJ,ktau);0;]);
            end
        end
        if stab
            Jc=Jc+prJc;
            Jctau=Jctau+prJctau;
        else
            Jc=Jc+prJc+prJctau;
        end
        %add partial matrices for LS phase condition
        if par.num.pctype~=0
            Mpresder1=Mpresder1+prMpresder1;
            Mpresder2=Mpresder2+prMpresder2;
        end
    end
    %trim matrix with the last raw being there only to keep sparse in the same partial size
    if par.num.pctype~=0
        Mpresder1=Mpresder1(1:end-1,1:end-1);
        Mpresder2=Mpresder2(1:end-1,1:end-1);
    end
    %boundary conditions
    if ~stab
        BIIR=reshape(permute(repmat(((1:NC2*P+1)-1)*(NC1*P+1)*n,[1 1 n])+repmat(permute(1:n,[1 3 2]),[1 NC2*P+1 1]),[1 3 2]),[(NC2*P+1)*n 1]);
        Jc((BIIR-1)*(Nom1*Nom2*n+2)+NC1*NC2*P^2*n+(1:length(BIIR)).')=1;
        BIIR=reshape(permute(repmat(((1:NC2*P+1)-1)*(NC1*P+1)*n+NC1*P*n,[1 1 n])+repmat(permute(1:n,[1 3 2]),[1 NC2*P+1 1]),[1 3 2]),[(NC2*P+1)*n 1]);
        Jc((BIIR-1)*(Nom1*Nom2*n+2)+NC1*NC2*P^2*n+(1:length(BIIR)).')=-1;
        BIIR=reshape(permute(repmat(((1:NC1*P)-1)*n,[1 1 n])+repmat(permute(1:n,[1 3 2]),[1 NC1*P 1]),[1 3 2]),[(NC1*P)*n 1]);
        Jc((BIIR-1)*(Nom1*Nom2*n+2)+NC1*NC2*P^2*n+(NC2*P+1)*n+(1:length(BIIR)).')=1;
        BIIR=reshape(permute((NC2*P)*(NC1*P+1)*n+repmat(((1:NC1*P)-1)*n,[1 1 n])+repmat(permute(1:n,[1 3 2]),[1 NC1*P 1]),[1 3 2]),[(NC1*P)*n 1]);
        Jc((BIIR-1)*(Nom1*Nom2*n+2)+NC1*NC2*P^2*n+(NC2*P+1)*n+(1:length(BIIR)).')=-1;
    end
end
if ~brute
    if ~stab
        %phase condition %not brute
        ukl=j2kl(pj(1:Nom1*Nom2*n));
        d1ukl=zeros(NC1*P+1,NC2*P+1,n);
        d1uklc=zeros(NC1,P+1,NC2*P+1,n);
        d2ukl=zeros(NC1*P+1,NC2*P+1,n);
        d2uklc=zeros(NC1,P+1,NC2*P+1,n);
        d1ucijc=zeros(NC1,P,NC2*P,n);
        d2ucijc=zeros(NC1,P,NC2*P,n);
        pfpt1c=zeros(NC1,P,NC2*P,n);
        pfpt2c=zeros(NC1,P,NC2*P,n);
        fd1ucijtauc=zeros(NC1,P,NC2*P,n);
        fd2ucijtauc=zeros(NC1,P,NC2*P,n);
        %PARFOR
        parfor kJ=1:NC1
            %             for kJ=1:NC1
            d1ukl_clJ=zeros(P+1,NC2*P+1,n);
            d2ukl_clJ=zeros(P+1,NC2*P+1,n);
            d1ucij_clJ=zeros(P,NC2*P,n);
            d2ucij_clJ=zeros(P,NC2*P,n);
            pfpt1_clJ=zeros(P,NC2*P,n);
            pfpt2_clJ=zeros(P,NC2*P,n);
            fd1ucijtau_clJ=zeros(P,NC2*P,n);
            fd2ucijtau_clJ=zeros(P,NC2*P,n);
            for kF=1:Ppc
                for lJ=1:NC2
                    pcijJ=zeros(Ppc,Ppc,npc);
                    pcijJtauom=zeros(Ppc,Ppc,npc,Ntau);
                    if kF==1
                        ukl_kl=repmat(permute(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),[4 5 1 2 3]),[P+1 P+1 1 1 1]);
                        %\partial u/\partial\theta_1 on t_{kl}
                        d1ukl_clJ(:,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td1Ppq,[1 1 1 1 n]),3),4));
                        %\partial u/\partial\theta_2 on t_{kl}
                        d2ukl_clJ(:,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td2Ppq,[1 1 1 1 n]),3),4));
                        %\partial u/\partial\theta_1 on c_{ij}
                        d1ucij_clJ(:,(lJ-1)*P+1:lJ*P,:)=cpstate_ex(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),numpars,1);
                        %\partial u/\partial\theta_2 on c_{ij}
                        d2ucij_clJ(:,(lJ-1)*P+1:lJ*P,:)=cpstate_ex(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),numpars,2);
                        %delayed state derivative w.r.t \theta_1
                        pcijJ=cpstate_ex(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),numpars,0);
                        for ktau=1:Ntau
                            pcijJtauom(:,:,:,ktau)=cptaustate_ex(ukl,om1*taus(ktau),om2*taus(ktau),(kJ-1)*Ppc+(1:(Ppc+1)),(lJ-1)*Ppc+(1:(Ppc+1)),numpars,0);
                        end
                    end
                    for lF=1:Ppc
                        %derivative w.r.t. explicite time related to omega 1
                        pfpt1_clJ(kF,(lJ-1)*P+lF,:)=par.system.sys_deri(c1ip((kJ-1)*Ppc+kF),c2jp((lJ-1)*Ppc+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1i);
                        %derivative w.r.t. explicite time related to omega 2
                        pfpt2_clJ(kF,(lJ-1)*P+lF,:)=par.system.sys_deri(c1ip((kJ-1)*Ppc+kF),c2jp((lJ-1)*Ppc+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,2i);
                        for ktau=1:Ntau
                            d1pcijJtauom=cptaustate_ex(ukl,om1*taus(ktau),om2*taus(ktau),(kJ-1)*Ppc+(1:(Ppc+1)),(lJ-1)*Ppc+(1:(Ppc+1)),numpars,1);
                            d2pcijJtauom=cptaustate_ex(ukl,om1*taus(ktau),om2*taus(ktau),(kJ-1)*Ppc+(1:(Ppc+1)),(lJ-1)*Ppc+(1:(Ppc+1)),numpars,2);
                            %derivative w.r.t. the delayed state
                            dfdutau=par.system.sys_deri(c1ip((kJ-1)*Ppc+kF),c2jp((lJ-1)*Ppc+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1+ktau);
                            %delayed vectorfield derivative w.r.t \theta_1
                            fd1ucijtau_clJ(kF,(lJ-1)*P+lF,:,ktau)=dfdutau*permute(d1pcijJtauom(kF,lF,:),[3 1 2]);
                            %delayed vectorfield derivative w.r.t \theta_2
                            fd2ucijtau_clJ(kF,(lJ-1)*P+lF,:,ktau)=dfdutau*permute(d2pcijJtauom(kF,lF,:),[3 1 2]);
                        end
                    end
                end
            end
            d1uklc(kJ,:,:,:)=d1ukl_clJ;
            d2uklc(kJ,:,:,:)=d2ukl_clJ;
            d1ucijc(kJ,:,:,:)=d1ucij_clJ;
            d2ucijc(kJ,:,:,:)=d2ucij_clJ;
            pfpt1c(kJ,:,:,:)=pfpt1_clJ;
            pfpt2c(kJ,:,:,:)=pfpt2_clJ;
            fd1ucijtauc(kJ,:,:,:)=fd1ucijtau_clJ;
            fd2ucijtauc(kJ,:,:,:)=fd2ucijtau_clJ;
        end
        %there is overlapping here so reshape is not possible
        for kJ=1:NC1
            d1ukl((kJ-1)*P+1:kJ*P+1,:,:)=d1ukl((kJ-1)*P+1:kJ*P+1,:,:)+permute(d1uklc(kJ,:,:,:),[2 3 4 1]);
            d2ukl((kJ-1)*P+1:kJ*P+1,:,:)=d2ukl((kJ-1)*P+1:kJ*P+1,:,:)+permute(d2uklc(kJ,:,:,:),[2 3 4 1]);
        end
        d1ucij=reshape(permute(d1ucijc,[2 1 3 4 5]),[NC1*P NC2*P n]);
        d2ucij=reshape(permute(d2ucijc,[2 1 3 4 5]),[NC1*P NC2*P n]);
        pfpt1=reshape(permute(pfpt1c,[2 1 3 4 5]),[NC1*P NC2*P n]);
        pfpt2=reshape(permute(pfpt2c,[2 1 3 4 5]),[NC1*P NC2*P n]);
        fd1ucijtau=reshape(permute(fd1ucijtauc,[2 1 3 4 5]),[NC1*P NC2*P n]);
        fd2ucijtau=reshape(permute(fd2ucijtauc,[2 1 3 4 5]),[NC1*P NC2*P n]);
        % phase conditions
        if par.num.pctype==0
            %Poincare-type phase condition <du,\Delta u>
            Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)+kl2j(d1ukl).';
            Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)+kl2j(d2ukl).';
        elseif par.num.pctype==1
            %uckl
            ucij=cstate(ukl);
            %d1Pcijkl
            Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)+(real(Mpresder1)*reshape(permute(repmat(wp.'*wp,[NC1 NC2 n]).*ucij,[1 2 3]),[NC1*P*NC2*P*n 1])).';
            Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)+(real(Mpresder2)*reshape(permute(repmat(wp.'*wp,[NC1 NC2 n]).*ucij,[1 2 3]),[NC1*P*NC2*P*n 1])).';
        end
        % frequency derivatives
        Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)+kl2j(d1ucij,NC1*P,NC2*P);
        Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)+kl2j(d2ucij,NC1*P,NC2*P);
        %explicite time dependency
        Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)+kl2j(pfpt1.*repmat(c1ip.',[1 NC2*P n])/om1^2,NC1*P,NC2*P);
        Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)+kl2j(pfpt2.*repmat(c2jp,[NC1*P 1 n])/om2^2,NC1*P,NC2*P);
        for ktau=1:Ntau
            Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+1)+kl2j(fd1ucijtau(:,:,:,ktau)*taus(ktau),NC1*P,NC2*P);
            Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)=Jc(1:NC1*NC2*P^2*n,Nom1*Nom2*n+2)+kl2j(fd2ucijtau(:,:,:,ktau)*taus(ktau),NC1*P,NC2*P);
        end
    end
else
    %phase condition %not brute
    ukl=j2kl(pj(1:Nom1*Nom2*n));
    %\partial u/\partial\theta_1 on t_{kl}
    d1ukl=zeros(NC1*P+1,NC2*P+1,n);
    %\partial u/\partial\theta_2 on t_{kl}
    d2ukl=zeros(NC1*P+1,NC2*P+1,n);
    for kJ=1:NC1
        for lJ=1:NC2
            ukl_kl=repmat(permute(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),[4 5 1 2 3]),[P+1 P+1 1 1 1]);
            %\partial u/\partial\theta_1 on t_{kl}
            d1ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td1Ppq,[1 1 1 1 n]),3),4));
            %\partial u/\partial\theta_2 on t_{kl}
            d2ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td2Ppq,[1 1 1 1 n]),3),4));
        end
    end
    Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+1,1:Nom1*Nom2*n)+kl2j(d1ukl).';
    Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)=Jc(Nom1*Nom2*n+2,1:Nom1*Nom2*n)+kl2j(d2ukl).';
end
% varargout

return;

    function uj=kl2j(varargin)
        if length(varargin)==1
            n1=Nom1;
            n2=Nom2;
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


    function pvcij=cpstate(pvij,deriv)
        %part approximation, size(pvkl)=[P+1 P+1 n], size(pcvkl)=[P P n]
        %declarate the result function
        %states at the collocation points
        %part of the representation points
        switch deriv
            case 0
                pvcij=permute(sum(sum(repmat(permute(pvij(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
                    .*repmat(cPpq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
            case 1
                pvcij=permute(sum(sum(repmat(permute(pvij(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
                    .*repmat(cd1Ppq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
            case 2
                pvcij=permute(sum(sum(repmat(permute(pvij(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
                    .*repmat(cd2Ppq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
        end
    end

    function vcij=cstate(vkl)
        %declarate the result function
        vcij=zeros(NC1*P,NC2*P,n);
        %states at the collocation points
        for ici=1:NC1
            for jci=1:NC2
                %part of the representation points
                vcij((ici-1)*P+1:ici*P,(jci-1)*P+1:jci*P,:)=cpstate(vkl((ici-1)*P+1:ici*P+1,(jci-1)*P+1:jci*P+1,:),0);
            end
        end
    end

end