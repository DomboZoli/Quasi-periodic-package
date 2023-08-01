function [cfoJac]=fun_cond_HKB_FFT(varargin)
% condition for torus locking based on relative harmonics determined by DFTs, not the best solution, evidently carries numerical issues, but good enough
%% INPUT
% pj=varargin{1};
% par=varargin{2};
% from here for jac_cond calc
% cpari=varargin{3};
% numpars=varargin{4};
% outG=varargin{5};
%% OUTPUT
%  cf: the condition
%  klmins: interval indeces in (k,l) for minima for all hi
%  minkl: iterated minima of composite chip thickness (2-gi)*hi*hi for all hi
%  th12min: iterated times for minima for all hi
%  himin: iterated actual minimum in hi
%  gimin: screen function at minima
%  phiimin: teeth position at the minima
%% VERSION v11.2(originated from fun_cond_HKB)
%   v1.0: simple NFO==1
%   v2.0: including the Jacobian too in cfoJac if numel(varargin)>2
%   v3.0: simple difference on the condition, make shure \omega_1>\omega_2
%   v4.0: |\omega_1-\omega_2|
%   v5.0: relative error in the spectrum
%   v6.0: complete velocity state
%   v7.0: vector field derivation taking into account the omega ukl dependency
%   v8.0: paralell
%   v8.2: minimum harmonics
%   v9.1: square difference
%   v9.1b: fft checking
%   v9.1c: smooth integration
%   v9.1d: takes vector field
%   v10.0: finished version
%   v11.0: paralellized
%   v11.1: TTwmin, TTwmax: included in nukmpars
%   v11.2: Dom^2
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% initial parameters
%initial state
pj=varargin{1};
par=varargin{2};
%allowable relative error on velocity spectra to the first harmonic
Dom=par.system.Dom;
%calculate the delay
taus=par.system.sys_tau([],par);
%number of delays
Ntau=numel(taus);
if ~isempty(varargin{3})
    %continuation parameters
    cpari=varargin{3};
    %Jacobian calculation
    Jaccalc=true;
else
    cpari=[];
    Jaccalc=false;
end
%numpars
numpars=varargin{4};
%wrinkle period T/Tw %
%min 
iTTwmin=numpars.condpar.iTTwmin;
%max 
iTTwmax=numpars.condpar.iTTwmax;
%data for external functions
Nom1=numpars.Nom1;
Nom2=numpars.Nom2;
n=numpars.n;
P=numpars.P;
NC1=numpars.NC1;
NC2=numpars.NC2;
%representation steps
h1=2*pi/(Nom1);
h2=2*pi/(Nom2);
%representation meshpoints
mth1=0:h1:(Nom1)*h1;
mth2=0:h2:(Nom2)*h2;
%shape functions
td1Ppq=numpars.td1Ppq;
td2Ppq=numpars.td2Ppq;
%period
T12=2*pi;
%frequency resolution
df12=1/T12;
%vibration frequencies
om1=real(pj(Nom1*Nom2*n+1,1));
om2=real(pj(Nom1*Nom2*n+2,1));
%define profile
ukl=j2kl(pj(1:Nom1*Nom2*n));
%determine phase condition
%\partial u/\partial\theta_1 on t_{kl}
d1ukl=zeros(NC1*P+1,NC2*P+1,n);
%\partial u/\partial\theta_2 on t_{kl}
d2ukl=zeros(NC1*P+1,NC2*P+1,n);
for i=1:NC1
    for j=1:NC2
        kJ=i;
        lJ=j;
        ukl_kl=repmat(permute(ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:),[4 5 1 2 3]),[P+1 P+1 1 1 1]);
        %\partial u/\partial\theta_1 on t_{kl}
        d1ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td1Ppq,[1 1 1 1 n]),3),4));
        %\partial u/\partial\theta_2 on t_{kl}
        d2ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:)=squeeze(sum(sum(ukl_kl.*repmat(td2Ppq,[1 1 1 1 n]),3),4));
    end
end
%ukltau
ukltau=zeros(NC1*P+1,NC2*P+1,n,Ntau);
% ucijtau=ctaustate(ukl,omT*tau,om*tau);
for ktau=1:Ntau
    ukltau(:,:,:,ktau)=taustate(ukl,om1*taus(ktau),om2*taus(ktau));
end
%determine the state of the milling system
fcij=par.system.sys_rhs(ukl,ukltau,ones(Nom1,Nom2),par,om1,om2);
%derivative by time
dukldt=fcij;%om1*d1ukl+om2*d2ukl;
%determine FFTs through \theta1 of velocities
dU1=fft(dukldt(1:end-1,:,:))/Nom1/df12;
%determine FFTs through \theta2 of velocities
dU2=fft(permute(dukldt(:,1:end-1,:),[2 1 3]))/Nom1/df12;
%first harmonic square to compare
S11=(reshape(dU1(2,1:Nom2-1,:),[(Nom2-1)*n 1])'*reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])+...
    reshape(dU1(2,2:Nom2,:),[(Nom2-1)*n 1])'*reshape(dU1(2,2:end,:),[(Nom2-1)*n 1]))/(4*pi);
S12=(reshape(dU2(2,1:Nom2-1,:),[(Nom2-1)*n 1])'*reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])+...
    reshape(dU2(2,2:Nom2,:),[(Nom2-1)*n 1])'*reshape(dU2(2,2:end,:),[(Nom2-1)*n 1]))/(4*pi);
%all harmonic from p+1 square to compare
Q1=(reshape(dU1(iTTwmin:iTTwmax,1:Nom2-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])+...
    reshape(dU1(iTTwmin:iTTwmax,2:Nom2,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))/(4*pi);
Q2=(reshape(dU2(iTTwmin:iTTwmax,1:Nom2-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])+...
    reshape(dU2(iTTwmin:iTTwmax,2:Nom2,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))/(4*pi);
%error
eps=S11+S12-Q1-Q2;
%relative error
epsilon=eps^2/(S11+S12)^2;
if ~Jaccalc
    %condition function for locking
    cfoJac=(epsilon-Dom^2);
else
    %check previous calculation
    P=numpars.P;
    cPpq=numpars.cPpq;
    tPpq=numpars.tPpq;
    %representation stepss
    h1=2*pi/(NC1*P);
    h2=2*pi/(NC2*P);
    %number of parameters
    np=numel(cpari);
    %number
    %% main code (no frequency fix)
    % representation interval length
    dth1=2*pi/NC1;
    dth2=2*pi/NC2;
    %derivative of ukl w.r.t. time
    duklkl=repmat(om1*td1Ppq(1:P,1:P,:,:)+om2*td2Ppq(1:P,1:P,:,:),[NC1 NC2 1 1]);
    %extend for periodicity
    duklkl(:,end+1,:,:)=duklkl(:,1,:,:);
    duklkl(end+1,:,:,:)=duklkl(1,:,:,:);
    %declare indeces
    i_1=zeros(1,NC1,Ntau);
    i_2=zeros(1,NC2,Ntau);
    Dc1=zeros(1,NC1,Ntau);
    Dc2=zeros(1,NC2,Ntau);
    %define delayed state on kl
    ukltau=zeros(Nom1,Nom2,n,Ntau);
    for ktau=1:Ntau
        for kJ=1:NC1
            %delayed mesh time (checked!)
            cijtau_1=mod((kJ-1)*dth1-om1*taus(ktau),2*pi);
            %calculate minimum modulo to a coll point
            i_1(:,kJ,ktau)=floor(cijtau_1/h1);
            %difference bitween the
            Dc1(:,kJ,ktau)=cijtau_1/h1-i_1(:,kJ,ktau);
        end
        for lJ=1:NC2
            %delayed mesh time (checked!)
            cijtau_2=mod((lJ-1)*dth2-om2*taus(ktau),2*pi);
            %calculate minimum modulo to a coll point
            i_2(:,lJ,ktau)=floor(cijtau_2/h2);
            %difference bitween the
            Dc2(:,lJ,ktau)=cijtau_2/h2-i_2(:,lJ,ktau);
        end
        ukltau(:,:,:,ktau)=taustate(ukl,om1*taus(ktau),om2*taus(ktau));
    end
    %uckltau
    ucijtau=zeros(NC1*P,NC2*P,n,Ntau);
    % ucijtau=ctaustate(ukl,omT*tau,om*tau);
    for ktau=1:Ntau
        ucijtau(:,:,:,ktau)=ctaustate(ukl,om1*taus(ktau),om2*taus(ktau));
    end
    %perform the relevant FFTs
    dU1klkl=fft(duklkl(1:end-1,:,:,:))/Nom1/df12;
    dU2klkl=fft(permute(duklkl(:,1:end-1,:,:),[2 1 3 4]))/Nom2/df12;
    dS11dukl=zeros(P+1,P+1);
    dS12dukl=zeros(P+1,P+1);
    dQ1dukl=zeros(P+1,P+1);
    dQ2dukl=zeros(P+1,P+1);
    for kP=1:P+1
        for lP=1:P+1
            %first harmonic square to compare
            dS11dukl(kP,lP)=2*(real(reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(repmat(dU1klkl(2,1:end-1,kP,lP),[1 1 n]),[(Nom2-1)*n 1]))+...
                real(reshape(dU1(2,2:end,:),[(Nom2-1)*n 1])'*reshape(repmat(dU1klkl(2,2:end,kP,lP),[1 1 n]),[(Nom2-1)*n 1])))/(4*pi);
            dS12dukl(kP,lP)=2*(real(reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(repmat(dU2klkl(2,1:end-1,kP,lP),[1 1 n]),[(Nom2-1)*n 1]))+...
                real(reshape(dU2(2,2:end,:),[(Nom2-1)*n 1])'*reshape(repmat(dU2klkl(2,2:end,kP,lP),[1 1 n]),[(Nom2-1)*n 1])))/(4*pi);
            %all harmonic from p+1 square to compare
            dQ1dukl(kP,lP)=2*(real(reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(repmat(dU1klkl(iTTwmin:iTTwmax,1:end-1,kP,lP),[1 1 n]),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
                real(reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(repmat(dU1klkl(iTTwmin:iTTwmax,2:end,kP,lP),[1 1 n]),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
            dQ2dukl(kP,lP)=2*(real(reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(repmat(dU2klkl(iTTwmin:iTTwmax,1:end-1,kP,lP),[1 1 n]),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
                real(reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(repmat(dU2klkl(iTTwmin:iTTwmax,2:end,kP,lP),[1 1 n]),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
        end
    end
    %determine Jacobian
    cfoJac=zeros(1,Nom1*Nom2*n+2+np);
    %this should be available
    dfdmuckJ=zeros(NC1,P+1,Nom2,n,np);
    %global Jac container, present state
    cfoJac_cpglo=zeros(NC1,NC2,(P+1)^2*n);
    %global Jac container, delayed state
    cfoJac_cdglo=zeros(NC1,NC2,(P+1)^2*n,Ntau);
    %     dfdmu=zeros(Nom1,Nom2,n,np);
    % TODO: paralell?
    parfor kJ=1:NC1
        %     for kJ=1:NC1
        %declaration of variables
        dS11dukl=zeros(P+1,P+1,n);
        dS12dukl=zeros(P+1,P+1,n);
        dQ1dukl=zeros(P+1,P+1,n);
        dQ2dukl=zeros(P+1,P+1,n);
        
        dS11dukltau=zeros(P+1,P+1,n,Ntau);
        dS12dukltau=zeros(P+1,P+1,n,Ntau);
        dQ1dukltau=zeros(P+1,P+1,n,Ntau);
        dQ2dukltau=zeros(P+1,P+1,n,Ntau);
        
        Gklkl=zeros(P+1,P+1,n,1,1,n);
        Gklkltau=zeros(P+1,P+1,n,1,1,n,Ntau);
        %define dfmu container
        dfdmucont=zeros(P+1,Nom2,n,np);
        %Jacobian container, present
        Jac_cpre=zeros(NC2,(P+1)^2*n);
        %Jacobian container, delayed
        Jac_cdel=zeros(NC2,(P+1)^2*n,Ntau);
        
        for lJ=1:NC2
            %                 if kF==1
            duklkl=zeros(Nom1,Nom2,n,P+1,P+1,n);
            duklkltau=zeros(Nom1,Nom2,n,P+1,P+1,n,Ntau);
            %cut out (k,l) (Ppc+1)^2*npc interval
            uklJ=ukl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:);
            %cut ukl tau
            uklJtau=ukltau((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:,:);
            %                 end
            
            for kF=1:P+1
                for lF=1:P+1
                    %derivative w.r.t. the state
                    Gklkl(kF,lF,:,1,1,:)=permute(par.system.sys_deri(mth1((kJ-1)*P+kF),mth2((lJ-1)*P+lF),permute(uklJ(kF,lF,:),[3,1,2]),permute(uklJtau(kF,lF,:,:),[3,4,2,1]),par,1),[3 4 1 5 6 2]);
                    for ktau=1:Ntau
                        Gklkltau(kF,lF,:,1,1,:,ktau)=permute(par.system.sys_deri(mth1((kJ-1)*P+kF),mth2((lJ-1)*P+lF),permute(uklJ(kF,lF,:),[3,1,2]),permute(uklJtau(kF,lF,:,:),[3,4,2,1]),par,1+ktau),[3 4 1 5 6 2 7]);
                    end
                    %derivatives w.r.t. the parameters
                    for knp=1:np
                        %                         dfdmu((kJ-1)*P+kF,(lJ-1)*P+lF,:,knp)=dfdmu((kJ-1)*P+kF,(lJ-1)*P+lF,:,knp)+permute(par.system.sys_deri(mth1((kJ-1)*P+kF),mth2((lJ-1)*P+lF),permute(uklJ(kF,lF,:),[3,1,2]),permute(uklJtau(kF,lF,:,:),[3,4,2,1]),par,2+cpari(knp)),[2 3 1]);
                        
                        dfdmucont(kF,(lJ-1)*P+lF,:,knp)=dfdmucont(kF,(lJ-1)*P+lF,:,knp)+permute(par.system.sys_deri(mth1((kJ-1)*P+kF),mth2((lJ-1)*P+lF),permute(uklJ(kF,lF,:),[3,1,2]),permute(uklJtau(kF,lF,:,:),[3,4,2,1]),par,2+cpari(knp)),[2 3 1]);
                    end
                end
            end
            %                 if kF==1
            %all needed due to FFT
            duklkl((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:,:,:,:)=repmat(Gklkl(:,:,:,1,1,:),[1 1 1 P+1 P+1 1]).*repmat(permute(tPpq(1:P+1,1:P+1,:,:),[1 2 5 3 4 6]),[1 1 n 1 1 n]);
            for ktau=1:Ntau
                duklkltau((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:,:,:,:,ktau)=duklkltau((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:,:,:,:,ktau)+repmat(Gklkltau(:,:,:,1,1,:,ktau),[1 1 1 P+1 P+1 1 1]).*repmat(permute(tPpq(1:P+1,1:P+1,:,:),[1 2 5 3 4 6 7]),[1 1 n 1 1 n 1]);
            end
            %cyclic boundary
            duklkl(:,end,:,:,:,:)=duklkl(:,1,:,:,:,:);
            duklkl(end,:,:,:,:,:)=duklkl(1,:,:,:,:,:);
            duklkltau(:,end,:,:,:,:,:)=duklkltau(:,1,:,:,:,:,:);
            duklkltau(end,:,:,:,:,:,:)=duklkltau(1,:,:,:,:,:,:);
            %perform the relevant FFTs (TODO: select important directions)
            dU1klkllJ=fft(duklkl(1:end-1,(lJ-1)*P+1:lJ*P+1,:,:,:,:))/Nom1/df12;
            dU2klklkL=fft(permute(duklkl((kJ-1)*P+1:kJ*P+1,1:end-1,:,:,:,:),[2 1 3 4 5 6]))/Nom1/df12;
            %perform the relevant FFTs
            dU1klkltaulJ=fft(duklkltau(1:end-1,(lJ-1)*P+1:lJ*P+1,:,:,:,:,:))/Nom1/df12;
            dU2klkltaukL=fft(permute(duklkltau((kJ-1)*P+1:kJ*P+1,1:end-1,:,:,:,:,:),[2 1 3 4 5 6 7]))/Nom1/df12;
            for k1=1:P+1
                for l1=1:P+1
                    for m1=1:n
                        %first harmonic square to compare
                        dS11dukl(k1,l1,m1)=2*(real(reshape(dU1(2,(kJ-1)*P+1:kJ*P,:),[(P)*n 1])'*reshape(dU1klkllJ(2,1:end-1,1:n,k1,l1,m1),[(P)*n 1]))+...
                            real(reshape(dU1(2,(kJ-1)*P+2:kJ*P+1,:),[(P)*n 1])'*reshape(dU1klkllJ(2,2:end,1:n,k1,l1,m1),[(P)*n 1])))/(4*pi);
                        dS12dukl(k1,l1,m1)=2*(real(reshape(dU2(2,(kJ-1)*P+1:kJ*P,:),[(P)*n 1])'*reshape(dU2klklkL(2,1:end-1,1:n,k1,l1,m1),[(P)*n 1]))+...
                            real(reshape(dU2(2,(kJ-1)*P+2:kJ*P+1,:),[(P)*n 1])'*reshape(dU2klklkL(2,2:end,1:n,k1,l1,m1),[(P)*n 1])))/(4*pi);
                        %all harmonic from p+1 square to compare
                        dQ1dukl(k1,l1,m1)=2*(real(reshape(dU1(iTTwmin:iTTwmax,(kJ-1)*P+1:kJ*P,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU1klkllJ(iTTwmin:iTTwmax,1:end-1,1:n,k1,l1,m1),[(iTTwmax-iTTwmin+1)*(P)*n 1]))+...
                            real(reshape(dU1(iTTwmin:iTTwmax,(kJ-1)*P+2:kJ*P+1,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU1klkllJ(iTTwmin:iTTwmax,2:end,1:n,k1,l1,m1),[(iTTwmax-iTTwmin+1)*(P)*n 1])))/(4*pi);
                        dQ2dukl(k1,l1,m1)=2*(real(reshape(dU2(iTTwmin:iTTwmax,(kJ-1)*P+1:kJ*P,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU2klklkL(iTTwmin:iTTwmax,1:end-1,1:n,k1,l1,m1),[(iTTwmax-iTTwmin+1)*(P)*n 1]))+...
                            real(reshape(dU2(iTTwmin:iTTwmax,(kJ-1)*P+2:kJ*P+1,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU2klklkL(iTTwmin:iTTwmax,2:end,1:n,k1,l1,m1),[(iTTwmax-iTTwmin+1)*(P)*n 1])))/(4*pi);
                        for ktau=1:Ntau
                            %first harmonic square to compare
                            dS11dukltau(k1,l1,m1,ktau)=2*(real(reshape(dU1(2,(kJ-1)*P+1:kJ*P,:),[(P)*n 1])'*reshape(dU1klkltaulJ(2,1:end-1,1:n,k1,l1,m1,ktau),[(P)*n 1]))+...
                                real(reshape(dU1(2,(kJ-1)*P+2:kJ*P+1,:),[(P)*n 1])'*reshape(dU1klkltaulJ(2,2:end,1:n,k1,l1,m1,ktau),[(P)*n 1])))/(4*pi);
                            dS12dukltau(k1,l1,m1,ktau)=2*(real(reshape(dU2(2,(kJ-1)*P+1:kJ*P,:),[(P)*n 1])'*reshape(dU2klkltaukL(2,1:end-1,1:n,k1,l1,m1,ktau),[(P)*n 1]))+...
                                real(reshape(dU2(2,(kJ-1)*P+2:kJ*P+1,:),[(P)*n 1])'*reshape(dU2klkltaukL(2,2:end,1:n,k1,l1,m1,ktau),[(P)*n 1])))/(4*pi);
                            %all harmonic from p+1 square to compare
                            dQ1dukltau(k1,l1,m1,ktau)=2*(real(reshape(dU1(iTTwmin:iTTwmax,(kJ-1)*P+1:kJ*P,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU1klkltaulJ(iTTwmin:iTTwmax,1:end-1,1:n,k1,l1,m1,ktau),[(iTTwmax-iTTwmin+1)*(P)*n 1]))+...
                                real(reshape(dU1(iTTwmin:iTTwmax,(kJ-1)*P+2:kJ*P+1,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU1klkltaulJ(iTTwmin:iTTwmax,2:end,1:n,k1,l1,m1,ktau),[(iTTwmax-iTTwmin+1)*(P)*n 1])))/(4*pi);
                            dQ2dukltau(k1,l1,m1,ktau)=2*(real(reshape(dU2(iTTwmin:iTTwmax,(kJ-1)*P+1:kJ*P,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU2klkltaukL(iTTwmin:iTTwmax,1:end-1,1:n,k1,l1,m1,ktau),[(iTTwmax-iTTwmin+1)*(P)*n 1]))+...
                                real(reshape(dU2(iTTwmin:iTTwmax,(kJ-1)*P+2:kJ*P+1,:),[(iTTwmax-iTTwmin+1)*(P)*n 1])'*reshape(dU2klkltaukL(iTTwmin:iTTwmax,2:end,1:n,k1,l1,m1,ktau),[(iTTwmax-iTTwmin+1)*(P)*n 1])))/(4*pi);
                        end
                    end
                end
            end
                                %set present Jacobian container
                                Jac_cpre(lJ,:)=reshape(2*eps/(S11+S12)^3*((dS11dukl+dS12dukl)*(Q1+Q2)-(dQ1dukl+dQ2dukl)*(S11+S12)),[(P+1)^2*n 1 1]).';
                                for ktau=1:Ntau
                                    %set delayed Jacobian container
                                    Jac_cdel(lJ,:,ktau)=reshape(2*eps/(S11+S12)^3*((dS11dukltau(:,:,:,ktau)+dS12dukltau(:,:,:,ktau))*(Q1+Q2)-(dQ1dukltau(:,:,:,ktau)+dQ2dukltau(:,:,:,ktau))*(S11+S12)),[(P+1)^2*n 1 1]).';
                                end
%             %present derivative
%             IRKJLJ=fun_IR(kJ,lJ,numpars);
%             %set Jacobian
%             cfoJac(1,IRKJLJ)=cfoJac(1,IRKJLJ)+reshape(2*eps/(S11+S12)^3*((dS11dukl+dS12dukl)*(Q1+Q2)-(dQ1dukl+dQ2dukl)*(S11+S12)),[(P+1)^2*n 1 1]).';
%             for ktau=1:Ntau
%                 %set weights
%                 inds1P2_w1=mod(i_1(1,kJ,ktau)+(0:P),NC1*P)+1;
%                 inds2P2_w1=mod(i_2(1,lJ,ktau)+(0:P),NC2*P)+1;
%                 %reindexed version
%                 DELEYIR=reshape(repmat((inds2P2_w1-1)*(NC1*P+1)*n,[P+1 1 n])+repmat(permute((inds1P2_w1-1)*n,[2 1 3]),[1 P+1 n])+repmat(permute(1:n,[1 3 2]),[P+1 P+1 1]),[(P+1)^2*n 1]);
%                 %set Jacobian
%                 cfoJac(1,DELEYIR)=cfoJac(1,DELEYIR)+reshape(2*eps/(S11+S12)^3*((dS11dukltau(:,:,:,ktau)+dS12dukltau(:,:,:,ktau))*(Q1+Q2)-(dQ1dukltau(:,:,:,ktau)+dQ2dukltau(:,:,:,ktau))*(S11+S12)),[(P+1)^2*n 1 1]).';
%             end
%                 end
        end
        dfdmuckJ(kJ,:,:,:,:)=dfdmucont;
        cfoJac_cpglo(kJ,:,:)=Jac_cpre;
        cfoJac_cdglo(kJ,:,:,:)=Jac_cdel;
    end

    
    dfdmu=zeros(Nom1,Nom2,n,np);
    %     for set correct position in the Jacobian
    for kJ=1:NC1
        %unfortunately there is overlapping and reshape cannot be used in dfmu
        dfdmu((kJ-1)*P+1:kJ*P+1,:,:,:)=dfdmu((kJ-1)*P+1:kJ*P+1,:,:,:)+permute(dfdmuckJ(kJ,:,:,:,:),[2 3 4 5 1]);
        for lJ=1:NC2
            %present indeces
            IRKJLJ=fun_IR(kJ,lJ,numpars);
            cfoJac(1,IRKJLJ)=cfoJac(1,IRKJLJ)+permute(cfoJac_cpglo(kJ,lJ,:),[2 3 1]);
            for ktau=1:Ntau
                %set weights
                inds1P2_w1=mod(i_1(1,kJ,ktau)+(0:P),NC1*P)+1;
                inds2P2_w1=mod(i_2(1,lJ,ktau)+(0:P),NC2*P)+1;
                %delayed indeces
                DELEYIR=reshape(repmat((inds2P2_w1-1)*(NC1*P+1)*n,[P+1 1 n])+repmat(permute((inds1P2_w1-1)*n,[2 1 3]),[1 P+1 n])+repmat(permute(1:n,[1 3 2]),[P+1 P+1 1]),[(P+1)^2*n 1]);
                cfoJac(1,DELEYIR)=cfoJac(1,DELEYIR)+permute(cfoJac_cdglo(kJ,lJ,:,ktau),[2 3 1 4]);
            end
        end
    end
    %set the periodicity on kl
    dfdmu(:,end,:,:)=dfdmu(:,1,:,:);
    dfdmu(end,:,:,:)=dfdmu(1,:,:,:);
    %velocity derivative w.r.t. frequency 1st coo
    dU1dom1=fft(d1ukl(1:end-1,:,:,:))/Nom1/df12;
    %velocity derivative w.r.t. frequency 2nd coo
    dU2dom1=fft(permute(d1ukl(:,1:end-1,:,:),[2 1 3 4]))/Nom1/df12;
    %first harmonic square to compare
    dS11dom1=2*(real(reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU1dom1(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU1(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU1dom1(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    dS12dom1=2*(real(reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU2dom1(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU2(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU2dom1(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    %all harmonic from p+1 square to compare
    dQ1dom1=2*(real(reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dom1(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dom1(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    dQ2dom1=2*(real(reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dom1(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dom1(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    %set jacobian
    cfoJac(1,Nom1*Nom2*n+1)=cfoJac(1,Nom1*Nom2*n+1)+2*eps/(S11+S12)^3*((Q1+Q2)*(dS11dom1+dS12dom1)-(dQ1dom1+dQ2dom1)*(S11+S12));
    %velocity derivative w.r.t. frequency 1st coo
    dU1dom2=fft(d2ukl(1:end-1,:,:,:))/Nom1/df12;
    %velocity derivative w.r.t. frequency 2nd coo
    dU2dom2=fft(permute(d2ukl(:,1:end-1,:,:),[2 1 3 4]))/Nom1/df12;
    %first harmonic square to compare
    dS11dom2=2*(real(reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU1dom2(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU1(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU1dom2(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    dS12dom2=2*(real(reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU2dom2(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU2(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU2dom2(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    %all harmonic from p+1 square to compare
    dQ1dom2=2*(real(reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dom2(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dom2(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    dQ2dom2=2*(real(reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dom2(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dom2(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    %set jacobian
    cfoJac(1,Nom1*Nom2*n+2)=cfoJac(1,Nom1*Nom2*n+2)+2*eps/(S11+S12)^3*((Q1+Q2)*(dS11dom2+dS12dom2)-(dQ1dom2+dQ2dom2)*(S11+S12));
    %derivatives w.r.t. the parameters
    %velocity derivative w.r.t. parameter 1 1st coo
    dU1dmu1=fft(dfdmu(1:end-1,:,:,1))/Nom1/df12;
    %velocity derivative w.r.t. parameter 1 2nd coo
    dU2dmu1=fft(permute(dfdmu(:,1:end-1,:,1),[2 1 3]))/Nom1/df12;
    %first harmonic square to compare
    dS11dmu1=2*(real(reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU1dmu1(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU1(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU1dmu1(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    dS12dmu1=2*(real(reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU2dmu1(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
        real(reshape(dU2(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU2dmu1(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
    %all harmonic from p+1 square to compare
    dQ1dmu1=2*(real(reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dmu1(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dmu1(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    dQ2dmu1=2*(real(reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dmu1(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
        real(reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dmu1(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
    %set jacobian
    cfoJac(1,Nom1*Nom2*n+2+1)=cfoJac(1,Nom1*Nom2*n+2+1)+2*eps/(S11+S12)^3*((Q1+Q2)*(dS11dmu1+dS12dmu1)-(dQ1dmu1+dQ2dmu1)*(S11+S12));
    if np>1
        %velocity derivative w.r.t. parameter 2 1st coo
        dU1dmu2=fft(dfdmu(1:end-1,:,:,2))/Nom1/df12;
        %velocity derivative w.r.t. parameter 2 2nd coo
        dU2dmu2=fft(permute(dfdmu(:,1:end-1,:,2),[2 1 3]))/Nom1/df12;
        %first harmonic square to compare
        dS11dmu2=2*(real(reshape(dU1(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU1dmu2(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
            real(reshape(dU1(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU1dmu2(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
        dS12dmu2=2*(real(reshape(dU2(2,1:end-1,:),[(Nom2-1)*n 1])'*reshape(dU2dmu2(2,1:end-1,1:n),[(Nom2-1)*n 1]))+...
            real(reshape(dU2(2,2:end,:),[(Nom2-1)*n 1])'*reshape(dU2dmu2(2,2:end,1:n),[(Nom2-1)*n 1])))/(4*pi);
        %all harmonic from p+1 square to compare
        dQ1dmu2=2*(real(reshape(dU1(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dmu2(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
            real(reshape(dU1(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU1dmu2(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
        dQ2dmu2=2*(real(reshape(dU2(iTTwmin:iTTwmax,1:end-1,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dmu2(iTTwmin:iTTwmax,1:end-1,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1]))+...
            real(reshape(dU2(iTTwmin:iTTwmax,2:end,:),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])'*reshape(dU2dmu2(iTTwmin:iTTwmax,2:end,1:n),[(iTTwmax-iTTwmin+1)*(Nom2-1)*n 1])))/(4*pi);
        %set jacobian
        cfoJac(1,Nom1*Nom2*n+2+2)=cfoJac(1,Nom1*Nom2*n+2+2)+2*eps/(S11+S12)^3*((Q1+Q2)*(dS11dmu2+dS12dmu2)-(dQ1dmu2+dQ2dmu2)*(S11+S12));
    end
    
end




return

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

    function pvcij=cpstate(pvij)
        %part approximation, size(pvkl)=[P+1 P+1 n], size(pcvkl)=[P P n]
        %declarate the result function
        %states at the collocation points
        %part of the representation points
        pvcij=permute(sum(sum(repmat(permute(pvij(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
            .*repmat(cPpq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
    end
%
    function pukl=j2kl(varargin)
        if length(varargin)==1
            pukl=permute(reshape(varargin{1},[n, Nom1, Nom2]),[2 3 1]);
        elseif length(varargin)==2
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{2}]),[2 3 1]);
        elseif length(varargin)==3
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{3}]),[2 3 1]);
        end
    end

    function vkltau=taustate(vkl,s1,s2)
        %extension part
        vkld=cat(2,vkl(1:Nom1-1,1:Nom2-1,:),vkl(1:Nom1-1,1:Nom2,:));
        vkldd=cat(1,vkld,vkld(1:Nom1-1,:,:));
        vkldd(2*Nom1-1,:,:)=vkldd(1,:,:);
        vkldd(:,2*Nom2-1,:)=vkldd(:,1,:);
        %         vkldd
        %the delay is produced in the representation points and interpolated at the original collocation points
        phis1=linspace(-2*pi,2*pi,Nom1*2-1);
        phis2=linspace(-2*pi,2*pi,Nom2*2-1);
        nr=size(vkl,3);
        vkltau=zeros(Nom1,Nom2,nr);
        for l_int=1:nr
            vkltau(:,:,l_int)=interp2(phis1,phis2,vkldd(:,:,l_int)',phis1(Nom1:end)-mod(s1,2*pi),(phis2(Nom2:end)-mod(s2,2*pi))','linear')';
        end
    end

end


