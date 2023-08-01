function eis = fun_invariance_2dtori(pj,par)
% build invariance equation for quasiperiod vibration
%% INPUT
%   pj: actual (jth) iteration step, [puij; pom;]
%       puij: (p: not iterated) representation points 
%       pom: (p: not iterated) vibration frequency
%   par: parameterlist (see parameterlist file)
%% OUTPUT
%   jc: colocation Jacobian
%% VERSION v4.0
%   v1.0: not optimized brute force solution to have a solution
%   v1.1: two unknown free frequencies
%   v1.1b: v1.1 modified for 1dof turning
%   v2.0: introduction of sys_rhs and sys_deri
%   v2.1: selectable scheme, par.num.scheme, 1: Chebishev, 2: Legendre
%   v3.0: idealistic case distinguished from flyover, par.system.NFO==0: ideal, par.system.NFO>0 flyover
%   v3.1: dim and dof in the par.system list, def_taus introduced
%   v3.5: possible to change between phase condition (pc) types: par.num.pctype==0: Poincaré pc, par.num.pctype==1: least square type
%   v4.0: finallized
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% preset parameters
%collocation order
P=par.num.p;
%collocation intervals for 2D tori both directions
NC1=par.num.NC1;
NC2=par.num.NC2;
%% main code
%set up phase condition reference solution
uj0=par.num.uj0;
%initialize integration points
cpps=zeros(7,7);
%initialize integration weights
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
%determine phase condition
%all shape functions in all collocation points
cPpq=zeros(P,P,P+1,P+1);
cd1Ppq=zeros(P,P,P+1,P+1);
cd2Ppq=zeros(P,P,P+1,P+1);
td1Ppq=zeros(P+1,P+1,P+1,P+1);
td2Ppq=zeros(P+1,P+1,P+1,P+1);
for pp=1:P+1
    for qq=1:P+1
        cPpq(:,:,pp,qq)=clag(pp,:)'*clag(qq,:);
        cd1Ppq(:,:,pp,qq)=cdlag(pp,:)'*clag(qq,:)/dth1;
        cd2Ppq(:,:,pp,qq)=clag(pp,:)'*cdlag(qq,:)/dth2;
        td1Ppq(:,:,pp,qq)=tdlag(pp,:)'*tlag(qq,:)/dth1;
        td2Ppq(:,:,pp,qq)=tlag(pp,:)'*tdlag(qq,:)/dth2;
    end
end

%% main code (milling related)
%calculate the delay
taus=par.system.sys_tau([],par);
%number of delays
Ntau=numel(taus);
%system DOF
dof=par.system.dof;
%dimension
n=2*dof;
%number of representation points
Nom1=NC1*P+1;
Nom2=NC2*P+1;
%interval points
th1i=((1:NC1+1)-1)*dth1;
% th2j=((1:NC2+1)-1)*dth2;
%collocation points
c1ip=reshape(repmat(th1i(1:end-1),[P,1])+repmat(cp'*dth1,[1 NC1]),[NC1*P,1])';
% c2jq=reshape(repmat(th2j(1:end-1),[P,1])+repmat(cp'*dth2,[1 NC2]),[NC2*P,1])';
%representation points
ukl=j2kl(pj(1:Nom1*Nom2*n));
%representation points uj0
ukl0=j2kl(uj0);
%vibration frequencies
om1=real(pj(Nom1*Nom2*n+1,1));
om2=real(pj(Nom1*Nom2*n+2,1));
%building of collocational equations
%(\patial u/\partial\theta_1)
d1ucij=zeros(NC1*P,NC2*P,n);
d2ucij=zeros(NC1*P,NC2*P,n);

d1ucij0=zeros(NC1*P,NC2*P,n);
d2ucij0=zeros(NC1*P,NC2*P,n);
for kJ=1:NC1
    for lJ=1:NC2
        ukl_ij=repmat(permute(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),[4 5 1 2 3]),[P P 1 1 1]);
        ukl_ij0=repmat(permute(ukl0((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),[4 5 1 2 3]),[P P 1 1 1]);
        
        d1ucij((kJ-1)*P+1:kJ*P,(lJ-1)*P+1:lJ*P,:)=sum(sum(ukl_ij.*repmat(cd1Ppq,[1 1 1 1 n]),3),4);
        d2ucij((kJ-1)*P+1:kJ*P,(lJ-1)*P+1:lJ*P,:)=sum(sum(ukl_ij.*repmat(cd2Ppq,[1 1 1 1 n]),3),4);
        
        d1ucij0((kJ-1)*P+1:kJ*P,(lJ-1)*P+1:lJ*P,:)=sum(sum(ukl_ij0.*repmat(cd1Ppq,[1 1 1 1 n]),3),4);
        d2ucij0((kJ-1)*P+1:kJ*P,(lJ-1)*P+1:lJ*P,:)=sum(sum(ukl_ij0.*repmat(cd2Ppq,[1 1 1 1 n]),3),4);
    end
end
%uckl
ucij=cstate(ukl);
%uckltau
ucijtau=zeros(NC1*P,NC2*P,n,Ntau);
% ucijtau=ctaustate(ukl,omT*tau,om*tau);
for ktau=1:Ntau
    ucijtau(:,:,:,ktau)=ctaustate(ukl,om1*taus(ktau),om2*taus(ktau));
end
% ucijtau=ctaustate(ukl,1*tau,om*tau);
%determine the state of the milling system
fcij=par.system.sys_rhs(ucij,ucijtau,repmat(c1ip.',[1 NC2*P]),par,om1,om2);
%dynamic part of the invariance equation
eis=kl2j(om1*d1ucij+om2*d2ucij-fcij,NC1*P,NC2*P);
% eis=kl2j(omT*d1ucij+om*d2ucij,NC1*P,NC2*P);
%boundary conditions
eis(NC1*P*NC2*P*n+1:NC1*P*NC2*P*n+Nom2*n)=kl2j(ukl(1,:,:)-ukl(end,:,:),1,Nom2);
eis(NC1*P*NC2*P*n+Nom2*n+1:NC1*P*NC2*P*n+Nom2*n+(Nom1-1)*n)=kl2j(ukl(1:end-1,1,:)-ukl(1:end-1,end,:),Nom1-1,1);
% par.num.pctype
if par.num.pctype==0
    %Poincare-type phase condition <du,\Delta u>
    %phase condition 1, intro
    eis(NC1*P*NC2*P*n+Nom2*n+(Nom1-1)*n+1)=0;
    %phase condition 2, intro
    eis(NC1*P*NC2*P*n+Nom2*n+(Nom1-1)*n+2)=0;
elseif par.num.pctype==1
    %least square phase condition with the reference solution being the current solution    
    %phase condition 1, intro
    eis(NC1*P*NC2*P*n+Nom2*n+(Nom1-1)*n+1)=reshape(d1ucij,[(NC1*P)^2*n 1])'*reshape(repmat(wp.'*wp,[NC1 NC2 n]).*ucij,[(NC1*P)^2*n 1]);
    %phase condition 2, intro
    eis(NC1*P*NC2*P*n+Nom2*n+(Nom1-1)*n+2)=reshape(d2ucij,[(NC1*P)^2*n 1])'*reshape(repmat(wp.'*wp,[NC1 NC2 n]).*ucij,[(NC1*P)^2*n 1]);
end

%end

return

    function pukl=j2kl(varargin)
        if length(varargin)==1
            pukl=permute(reshape(varargin{1},[n, Nom1, Nom2]),[2 3 1]);
        elseif length(varargin)==2
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{2}]),[2 3 1]);
        elseif length(varargin)==3
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{3}]),[2 3 1]);
        end
    end

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


    function vcijtau=ctaustate(vkl,s1,s2)
        %extension part
        vkld=cat(2,vkl(1:Nom1-1,1:Nom2-1,:),vkl(1:Nom1-1,1:Nom2,:));
        vkldd=cat(1,vkld,vkld(1:Nom1-1,:,:));
        vkldd(2*Nom1-1,:,:)=vkldd(1,:,:);
        vkldd(:,2*Nom2-1,:)=vkldd(:,1,:);
%         vkldd
        %the delay is produced in the representation points and interpolated at the original collocation points 
        phis1=linspace(-2*pi,2*pi,Nom1*2-1);
        phis2=linspace(-2*pi,2*pi,Nom2*2-1);
        vkltau=zeros(Nom1,Nom2,n);
        for l_int=1:n
                vkltau(:,:,l_int)=interp2(phis1,phis2,vkldd(:,:,l_int)',phis1(Nom1:end)-mod(s1,2*pi),(phis2(Nom2:end)-mod(s2,2*pi))','linear')';
        end
        %declarate the result function
        vcijtau=zeros(NC1*P,NC2*P,n);
        %states at the collocation points
        for ici=1:NC1
            for jci=1:NC2
                %part of the representation points
                vcijtau((ici-1)*P+1:ici*P,(jci-1)*P+1:jci*P,:)=cpstate(vkltau((ici-1)*P+1:ici*P+1,(jci-1)*P+1:jci*P+1,:));
            end
        end
    end

    function vcij=cstate(vkl)
        %declarate the result function
        vcij=zeros(NC1*P,NC2*P,n);
        %states at the collocation points
        for ici=1:NC1
            for jci=1:NC2
                %part of the representation points
                vcij((ici-1)*P+1:ici*P,(jci-1)*P+1:jci*P,:)=cpstate(vkl((ici-1)*P+1:ici*P+1,(jci-1)*P+1:jci*P+1,:));
            end
        end
    end

    function pvcij=cpstate(pvij)
        %part approximation, size(pvkl)=[P+1 P+1 n], size(pcvkl)=[P P n]
        %declarate the result function
        %states at the collocation points
        %part of the representation points
        pvcij=permute(sum(sum(repmat(permute(pvij(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
            .*repmat(cPpq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
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


end

