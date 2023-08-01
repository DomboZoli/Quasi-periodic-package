function Jres=jac_cpari_ex(pj,par,cpari,outG,numpars)
%function to provide the derivative of invariant tori equation w.r.t. the
%parameters, it might needed to conditions
%% VERSION 2.0
%v1.0: base code
%v2.0: parallelized
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
% main code: outsourced due to paralellization
%data for external functions
Nom1=numpars.Nom1;
Nom2=numpars.Nom2;
P=numpars.P;
n=numpars.n;
cPpq=numpars.cPpq;
cd1Ppq=numpars.cd1Ppq;
cd2Ppq=numpars.cd2Ppq;
np=numel(cpari);
NC1=numpars.NC1;
NC2=numpars.NC2;
c1ip=numpars.c1ip;
c2jp=numpars.c2jp;

%preallocate Jacobian structure
Jres=spalloc(Nom1*Nom2*n+2+np,np,ceil(0.1*(Nom1*Nom2*n+2)^2));
%building the Jacobian w.r.t. the parameters
ipar=par;
%frequencies
om1=pj(Nom1*Nom2*n+1);
om2=pj(Nom1*Nom2*n+2);
%derive pkl
pkl_cp=j2kl(pj(1:Nom1*Nom2*n),Nom1,Nom2,n);
%calculate the delays
[taupcs,npt]=par.system.sys_tau([],par);
%number of delays
Ntau=numel(taupcs);
%delay indeces
%         [~,npt]=par.system.sys_tau([],par);
%prepare derivatives
for knp=1:np
    %delay parameter index
    dpi=find(npt==cpari(knp),1);
    %the parameter is a delay?
    knpisdelay=~isempty(dpi);
    %             knpisdelay=false;
    %check if precalculate Gtau is available
    if knpisdelay && ~isempty(outG)
        outGvar=true;
        outGtau=outG{1};
    else
        outGvar=false;
        outGtau=[];
    end
     dmuknpc=zeros(NC1,NC2,P,P,n);
%     tic;
    %PARFOR
    parfor kJ=1:NC1
%         for kJ=1:NC1
                    %allocate the derivative
            dmuknp=zeros(NC2,P,P,n);
        for lJ=1:NC1
            %cut out (k,l) (Ppc+1)^2*npc interval
            pklJ=pkl_cp((kJ-1)*P+1:kJ*P+1,(lJ-1)*P+1:lJ*P+1,:);
            %present state
            pcijJ=cpstate_ex(pklJ,numpars,0);
            %delayed state
            pcijJtauom=zeros(P,P,n,Ntau);
            %delayed state
            for ktau=1:Ntau
                pcijJtauom(:,:,:,ktau)=cptaustate_ex(pkl_cp,om1*taupcs(ktau),om2*taupcs(ktau),(kJ-1)*P+(1:(P+1)),(lJ-1)*P+(1:(P+1)),numpars,0);
            end
            %force derivative w.r.t xtau
            Gijkltau=zeros(P,P,n,1,1,n,Ntau);
            if knpisdelay
                if ~outGvar
                    for kF=1:P
                        for lF=1:P
                            for ktau=1:Ntau
                                %TODO: for multiple delay case only the
                                %dpi should be derived, should be
                                %pulled from the main code
                                Gijkltau(kF,lF,:,1,1,:,ktau)=permute(par.system.sys_deri(c1ip((kJ-1)*P+kF),c2jp((lJ-1)*P+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1+ktau),[3,4,1,5,6,2]);
                            end
                        end
                    end
                else
                    Gijkltau(:,:,:,1,1,:,:)=outGtau(:,:,:,kJ,lJ,:,1:Ntau);
                end
                %delayed derivative state
                for ktau=1:Ntau
                    %normal derivative w.r.t. the delay parameter
                    pcijJtauom(:,:,:,ktau)=om1*cptaustate_ex(pkl_cp,om1*taupcs(ktau),om2*taupcs(ktau),(kJ-1)*P+(1:(P+1)),(lJ-1)*P+(1:(P+1)),numpars,1)+...
                        om2*cptaustate_ex(pkl_cp,om1*taupcs(ktau),om2*taupcs(ktau),(kJ-1)*P+(1:(P+1)),(lJ-1)*P+(1:(P+1)),numpars,2);
                end
            end

            for kF=1:P
                for lF=1:P
                    %TODO: check this part, continuation was not tried for multple delays
                    %this is not general for general delays
                    if ~knpisdelay
                        %the knp parameter is not a delay
                        dmuknp(lJ,kF,lF,:)=-par.system.sys_deri(c1ip((kJ-1)*P+kF),c2jp((lJ-1)*P+lF),permute(pcijJ(kF,lF,:),[3,1,2]),permute(pcijJtauom(kF,lF,:,:),[3,4,2,1]),par,1+Ntau+cpari(knp));
                    else
                        %TODO: make it efficient, the knp parameter is a delay (delayed)
                        %                                 dmuknp(kF,lF,:)=squeeze(sum(sum(sum(repmat(Gijkltau(kF,lF,:,1,1,:,dpi),[1 1 1 P P 1 1]).*repmat(permute(pcijJtauom(:,:,:,dpi),[4 5 6 1 2 3]),[1 1 n 1 1 1]),6),4),5));
                        dmuknp(lJ,kF,lF,:)=dmuknp(lJ,kF,lF,:)+permute(permute(Gijkltau(kF,lF,:,1,1,:,dpi),[3 6 1 2 4 5 7])*squeeze(pcijJtauom(kF,lF,:,dpi)),[2 3 4 1]);
                    end
                end
            end

        end
        dmuknpc(kJ,:,:,:,:)=dmuknp;
    end
    %storing data
    for kJ=1:NC1
        for lJ=1:NC2
                %                     col_der_np(:,kJ,lJ)=reshape(dmuknp,[P^2*n 1]);
            %store in jacobian
            spi=zeros(P^2*n+1,1);
            spi(P^2*n+1,1)=Nom1*Nom2*n+2+np;
            %collocation indeces
            spi(1:P^2*n,1)=fun_IC(kJ,lJ,numpars);
            %Jacobian
            Jres(1:Nom1*Nom2*n+2+np,knp)=Jres(1:Nom1*Nom2*n+2+np,knp)+sparse(spi,1,[reshape(squeeze(dmuknpc(kJ,lJ,:,:,:)),[P^2*n 1]);0;]);
        end
    end
%     toc
end
return