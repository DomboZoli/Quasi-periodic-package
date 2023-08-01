function handle=br_plot(branch,clr,type)
%v2.1: reject empty points
%v3.0: come back with plot handle
%v4.0: imag(): dimension
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code
%colors
colors='brmgy';
%BR_PLOT Summary of this function goes here
% fit for collocation method
%   Detailed explanation goes here
if real(type)==1
    %only norm
    pbl=zeros(length(branch.points),2);
    l=1;
    for k=1:length(branch.points)
        if ~isempty(branch.points(k).omega1)
            pbl(l,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(l,2)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
            l=l+1;
        end
    end
    handle=plot(pbl(1:l-1,1),pbl(1:l-1,2),clr);
elseif real(type)==1.1
    %peak 2 peak amplitude
    %only norm
    pbl=zeros(length(branch.points),2);
    l=1;
    for k=1:length(branch.points)
        if ~isempty(branch.points(k).omega1)
            pbl(l,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(l,2)=max(max(branch.points(k).profile(:,:,1),[],1),[],2)-min(min(branch.points(k).profile(:,:,1),[],1),[],2);
            l=l+1;
        end
    end
    handle=plot(pbl(1:l-1,1),pbl(1:l-1,2),clr);
elseif real(type)==1.2
    %based on expected rotation number automatically given by the omega1/omega2
    %only rough discrete picking form profile
    pbl=zeros(length(branch.points),2);
    l=1;
    Nx=size(branch.points(1).profile,1);
    Ny=size(branch.points(1).profile,2);
    for k=1:numel(branch.points)
        if ~isempty(branch.points(k).omega1)
            %rotation number
            rl=branch.points(k).omega1/branch.points(k).omega2;
            
            xi=mod(rl*linspace(1,Nx,Nx)-1,Nx)+1;
            yi=mod((1/rl)*linspace(1,Ny,Ny)-1,Ny)+1;
            ri=round((yi-1)*Ny+xi);
            pr1=branch.points(k).profile(:,:,1);
            pbl(l,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(l,2)=max(pr1(ri),[],2)-min(pr1(ri),[],2);
            l=l+1;
        end
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
elseif real(type)==1.3
    %based on expected rotation number automatically given by the omega1/omega2
    %2D interpolated
    pbl=zeros(length(branch.points),2);
    l=1;
    for k=1:numel(branch.points)
        if ~isempty(branch.points(k).omega1)
            %rotation number
            rl=branch.points(k).omega1/branch.points(k).omega2;
            xi=mod(rl*branch.mesh.m1,2*pi);
            yi=mod((1/rl)*branch.mesh.m2,2*pi);
            pr1ri=interp2(branch.mesh.m1,branch.mesh.m2,branch.points(k).profile(:,:,1)',xi,yi,'spline')';
            pbl(l,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(l,2)=max(pr1ri)-min(pr1ri);
            l=l+1;
        end
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
elseif real(type)==2
    %Soboljev norm
    pbl=zeros(length(branch.points),2);
    for k=1:length(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=norm(kl2j(branch.points(k).profile(:,:,:),size(branch.points(k).profile,1)),2);
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
elseif real(type)==3
    %vibration frequency 1
    pbl=zeros(length(branch.points),2);
    for k=1:length(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=branch.points(k).omega1;
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
elseif real(type)==3.5
    %vibration frequency 2
    pbl=zeros(length(branch.points),2);
    for k=1:length(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=branch.points(k).omega2;
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
    
elseif real(type)==4
    %only norm with index numbers
    pbl=zeros(2,length(branch.points));
    l=1;
    for k=1:length(branch.points)
        if ~isempty(branch.points(k).omega1)
            l=l+1;
            if length(branch.contpar)==1
                pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
                pbl(k,2)=max(max(branch.points(k).profile(:,:,1),[],1),[],2)-min(min(branch.points(k).profile(:,:,1),[],1),[],2);
            else
                pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(1)});
                pbl(k,2)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(2)});
            end
        end
    end
    handle=plot(pbl(1:l-1,1),pbl(1:l-1,2),clr);
    hold on;
    for k=1:l-1
        text(pbl(k,1),pbl(k,2),['---' int2str(k)]);
    end
elseif real(type)==5
    %plot rotation number
    pbl=zeros(2,length(branch.points));
    rot=zeros(1,length(branch.points));
    for k=1:length(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
        rot(k,1)=branch.points(k).omega1/(branch.points(k).omega2);
        
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
    hold on;
    for k=1:3:length(branch.points)
        text(pbl(k,1),pbl(k,2),['-------' num2str(rot(k))]);
    end
elseif real(type)==5.1
    %plot rotation number, no text
    pbl=zeros(2,numel(branch.points));
    rot=zeros(1,numel(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
        rot(k,1)=branch.points(k).omega1/(branch.points(k).omega2);
        
    end
    handle=plot(pbl(:,1),rot(:,1),clr);
    %         hold on;
    %         for k=1:3:length(branch.points)
    %             text(pbl(k,1),pbl(k,2),['-------' num2str(rot(k))]);
    %         end
elseif real(type)==5.2
    %plot rotation number, no text
    pbl=zeros(2,numel(branch.points));
    rot=zeros(1,numel(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
        pbl(k,2)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
        rot(k,1)=branch.points(k).omega2/(branch.points(k).omega1);
        
    end
    handle=plot(pbl(:,1),rot(:,1),clr);
    %         hold on;
    %         for k=1:3:length(branch.points)
    %             text(pbl(k,1),pbl(k,2),['-------' num2str(rot(k))]);
    %         end
elseif real(type)==6
    %two params
    pbl=zeros(2,length(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(2)});
        pbl(k,2)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(1)});
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
elseif real(type)==6.1
    %two params (reversed)
    pbl=zeros(2,length(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(2)});
        pbl(k,2)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(1)});
    end
    handle=plot(pbl(:,2),pbl(:,1),clr);
elseif real(type)==6.5
    %two params
    pbl=zeros(2,length(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(2)});
        pbl(k,2)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(1)});
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
    hold on;
    for k=1:numel(branch.points)
        text(pbl(k,1),pbl(k,2),['-------' num2str(k)]);
    end
elseif real(type)==6.6
    %two params
    pbl=zeros(2,length(branch.points));
    for k=1:numel(branch.points)
        pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(2)});
        pbl(k,2)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar(1)});
    end
    handle=plot(pbl(:,2),pbl(:,1),clr);
    hold on;
    for k=1:numel(branch.points)
        text(pbl(k,2),pbl(k,1),['-------' num2str(k)]);
    end
elseif real(type)==10
    %only norm with index numbers
    pbl=zeros(2,length(branch.points));
    l=1;
    for k=1:length(branch.points)
        if ~isempty(branch.points(k)) && ~isempty(branch.points(k).profile)
            pbl(1,l)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(2,l)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
            l=l+1;
        end
    end
    Np=l-1;
    handle=plot(pbl(1,:),pbl(2,:),clr);
    ll=1;
    for k=1:Np
        if ~isempty(branch.points(k)) && ~isempty(branch.points(k).profile)
            par=branch.points(k).par;
            % collocation intervals
            NC1=par.num.NC1;
            NC2=par.num.NC2;
            %polynomial order
            P=par.num.p;
            % Nom=10;
            Nom1=NC1*P+1;
            % discretization parameter for the vibration frequency
            Nom2=NC2*P+1;
            %representation steps
            h1=2*pi/(NC1*P);
            h2=2*pi/(NC2*P);
            % tooth passing frequency
            omT=branch.points(k).par.system.Om*branch.points(k).par.system.Z;
            Om=branch.points(k).par.system.Om;
            %size
            nm=length(branch.points(k).par.system.xis);
            %delay
            tau=2*pi/omT;
            % vibration frequency
            om=branch.points(k).omega2;
            %profile
            v0kl=branch.points(k).profile(:,:,1:nm);
            v0kld=cat(2,v0kl,v0kl(:,2:Nom2,:));
            v0kldd=cat(1,v0kld,v0kld(2:Nom1,:,:));
            %span
            phis1=linspace(-2*pi,2*pi,Nom1*2-1);
            phis2=linspace(-2*pi,2*pi,Nom2*2-1);
            %preallocation
            v0tau=zeros(Nom1,Nom2,nm);
            %interpolation
            for l_int=1:nm
                v0tau(:,:,l_int)=interp2(phis1,phis2,v0kldd(:,:,l_int)',phis1(Nom1:end)-mod(omT*tau,2*pi),(phis2(Nom2:end)-mod(om*tau,2*pi))','nearest')';
            end
            Dx=zeros(Nom1,Nom2,3);
            for ln=1:nm
                Dx(:,:,1)=Dx(:,:,1)+par.system.us(1,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
                Dx(:,:,2)=Dx(:,:,2)+par.system.us(2,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
                Dx(:,:,3)=Dx(:,:,3)+par.system.us(3,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
            end
            %representation meshpoints
            mth1=0:h1:(NC1*P)*h1;
            mth2=0:h2:(NC2*P)*h2;
            %leaves cut
            lc=1;
            for l=1:branch.points(k).par.system.Z
                phil=repmat(Om/omT*mth1.',[1 NC2*P+1])+(l-1)*2*pi/branch.points(k).par.system.Z;
                htl=double(0<((branch.points(k).par.system.fZ+Dx(:,:,1)).*sin(phil).*sin(par.system.kappa)+Dx(:,:,2).*cos(phil)*sin(par.system.kappa)))-Dx(:,:,3)*cos(par.system.kappa);
                %                 gi0_f=double((branch.points(k).par.system.phien+dphi)...
                %                     <=mod(phil,2*pi) & mod(phil,2*pi)<=...
                %                     branch.points(k).par.system.phiex-dphi);
                %                     dphi=branch.points(k).par.system.Dphi;
                dphi=0*2.43/180*pi;
                gi0_f=double((branch.points(k).par.system.phien+dphi)...
                    <=mod(phil,2*pi) & mod(phil,2*pi)<=...
                    branch.points(k).par.system.phiex-dphi);
                shl=prod(prod(htl(gi0_f>0),1),2);
                lc=lc*shl;
            end
            if lc==0
                %one edge left the surface
                hold on;
                plot(branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar}),max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2),'ro');
            end
            ll=ll+1;
        end
    end
elseif real(type)==11
    %condition function
    pbl=zeros(length(branch.points),2);
    for k=1:length(branch.points)
        if ~isempty(branch.points(k).omega1)
            pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            pbl(k,2)=max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2);
        end
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
    for k=1:length(branch.points)
        uj=kl2j(branch.points(k).profile,branch.points(k).par.num.NC1*branch.points(k).par.num.p+1);
        uj(end+1,1)=branch.points(k).omega1;
        uj(end+1,1)=branch.points(k).omega2;
        [cf,~,~]=fun_condition(uj,branch.points(k).par);
        if cf>0
            hold on;
            plot(branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar}),max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2),'ko');
        else
            hold on;
            plot(branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar}),max(max(abs(branch.points(k).profile(:,:,1)),[],1),[],2),'mo');
        end
    end
    
elseif real(type)==12
    %plot Ast+A diagram
    pbl=zeros(3,length(branch.points));
    l=1;
    par_def=def_Z2_NCT_Alb_l_2dof_dphi_c();
    %         par_def=def_SV_selmodes_P10_xyz_b_real_en_modal();
    %         par_def=def_SV_selmodes_P9_xyz_b_m_real4_18_En_modal();
    par_def.num.parameters(2)=1;
    for k=1:length(branch.points)
        if ~isempty(branch.points(k)) && ~isempty(branch.points(k).profile)
            pbl(1,l)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            %calculate the stationary solution, WARNING: only for specific def files
            par_def.process.parameters=[branch.points(k).par.system.ap 2*pi/branch.points(k).par.system.Om];
            par_def.cutter.parameters(4)=branch.points(k).par.system.ap;
            par_def=zd_geomserrcut(par_def);
            N1=size(branch.points(k).profile(:,:,1),1)-1;
            %calculate the stationary solution
            stat_sol=zd_statlinsol_deffile(par_def,N1);
            
            %stationary amplitude, x direction
            %                 pbl(3,l)=max(stat_sol.xst(:,1))-min(stat_sol.xst(:,1))/2;
            pbl(3,l)=2*sqrt(stat_sol.XST(2,2)*stat_sol.XST(2,2)')*1e6;
            %determine the profile regarding to the second frequency
            x2=zeros([size(branch.points(k).profile,1) size(branch.points(k).profile,2) 3]);
            X2=zeros([size(branch.points(k).profile,1)-1 size(branch.points(k).profile,2) 3]);
            
            for sk=1:3
                for lk=1:size(branch.points(k).par.system.us,2)
                    x2(:,:,sk)=x2(:,:,sk)+branch.points(k).par.system.us(sk,lk)*(branch.points(k).profile(:,:,lk));
                end
                
            end
            %remove the stationary part
            %                     x2(:,:,sk)=x2(:,:,sk)-(repmat(real([stat_sol.qst(:,lk);stat_sol.qst(1,lk);]),[1 size(branch.points(k).profile,2)])));
            x2=x2-repmat(permute(interp1([stat_sol.tm.';stat_sol.period;]/stat_sol.period,[stat_sol.xst;stat_sol.xst(1,:);],linspace(0,1,size(branch.points(k).profile,1)).','linear','extrap'),[1 3 2]),[1 size(branch.points(k).profile,2) 1]);
            %theoretically it can be only real
            x2=real(x2);
            for sk=1:3
                X2(:,:,sk)=fft(x2(:,1:end-1,sk).',size(x2,1)-1)/(size(x2,1)-1);
            end
            X2=permute(sum(X2,2)/size(x2,1),[1 3 2]);
            %the amplitude
            %average
            pbl(2,l)=2*sqrt(X2(2,2)*X2(2,2)')*1e6;
            l=l+1;
        end
    end
    Np=l-1;
    %stationary solution
    hold on
    handle=plot(pbl(1,:)*1000,abs(pbl(3,:)),clr);
    %Ast_1+A_1
    plot(pbl(1,:)*1000,abs(pbl(2,:))+abs(pbl(3,:)),clr);
    ll=1;
    
    for k=1:Np
        if ~isempty(branch.points(k)) && ~isempty(branch.points(k).profile)
            ll=ll+1;
            par=branch.points(k).par;
            % collocation intervals
            NC1=par.num.NC1;
            NC2=par.num.NC2;
            %polynomial order
            P=par.num.p;
            % Nom=10;
            Nom1=NC1*P+1;
            % discretization parameter for the vibration frequency
            Nom2=NC2*P+1;
            %representation steps
            h1=2*pi/(NC1*P);
            h2=2*pi/(NC2*P);
            % tooth passing frequency
            omT=branch.points(k).par.system.Om*branch.points(k).par.system.Z;
            Om=branch.points(k).par.system.Om;
            %size
            nm=length(branch.points(k).par.system.xis);
            %delay
            tau=2*pi/omT;
            % vibration frequency
            om=branch.points(k).omega2;
            %profile
            v0kl=branch.points(k).profile(:,:,1:nm);
            v0kld=cat(2,v0kl,v0kl(:,2:Nom2,:));
            v0kldd=cat(1,v0kld,v0kld(2:Nom1,:,:));
            %span
            phis1=linspace(-2*pi,2*pi,Nom1*2-1);
            phis2=linspace(-2*pi,2*pi,Nom2*2-1);
            %preallocation
            v0tau=zeros(Nom1,Nom2,nm);
            %interpolation
            for l_int=1:nm
                v0tau(:,:,l_int)=interp2(phis1,phis2,v0kldd(:,:,l_int)',phis1(Nom1:end)-mod(omT*tau,2*pi),(phis2(Nom2:end)-mod(om*tau,2*pi))','nearest')';
            end
            Dx=zeros(Nom1,Nom2,3);
            for ln=1:nm
                Dx(:,:,1)=Dx(:,:,1)+par.system.us(1,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
                Dx(:,:,2)=Dx(:,:,2)+par.system.us(2,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
                Dx(:,:,3)=Dx(:,:,3)+par.system.us(3,ln)*(v0kl(:,:,ln)-v0tau(:,:,ln));
            end
            %representation meshpoints
            mth1=0:h1:(NC1*P)*h1;
            mth2=0:h2:(NC2*P)*h2;
            %leaves cut
            lc=1;
            for l=1:branch.points(k).par.system.Z
                phil=repmat(Om/omT*mth1.',[1 NC2*P+1])+(l-1)*2*pi/branch.points(k).par.system.Z;
                htl=double(0<((branch.points(k).par.system.fZ+Dx(:,:,1)).*sin(phil).*sin(par.system.kappa)+Dx(:,:,2).*cos(phil)*sin(par.system.kappa)))-Dx(:,:,3)*cos(par.system.kappa);
                %                 gi0_f=double((branch.points(k).par.system.phien+dphi)...
                %                     <=mod(phil,2*pi) & mod(phil,2*pi)<=...
                %                     branch.points(k).par.system.phiex-dphi);
                %                     dphi=branch.points(k).par.system.Dphi;
                dphi=2.5/180*pi;
                %                     dphi=10/180*pi;
                gi0_f=double(max(branch.points(k).par.system.phien,dphi)...
                    <=mod(phil,2*pi) & mod(phil,2*pi)<=...
                    min(branch.points(k).par.system.phiex,pi-dphi));
                shl=prod(prod(htl(gi0_f>0),1),2);
                lc=lc*shl;
            end
            if lc==0
                %one edge left the surface
                hold on;
                plot(pbl(1,ll)*1000,abs(pbl(2,ll))+abs(pbl(3,ll)),'ro');
                %                      text(pbl(1,ll)*1000,abs(pbl(2,ll))+abs(pbl(3,ll)),['---' int2str(k)]);
            end
        end
    end
elseif real(type)==99
    %scaled norm
    pbl=zeros(length(branch.points),2);
    for k=1:length(branch.points)
        if ~isempty(branch.points(k).omega1)
            pbl(k,1)=branch.points(k).par.system.(branch.points(k).par.parlist{branch.contpar});
            point = branch.points(k);
            tau=2*pi/point.par.system.Om;
            deltax=max(max(point.profile(:,:,1),[],1),[],2)-...
                min(min(point.profile(:,:,1),[],1),[],2);
            deltaxdot=max(max(point.profile(:,:,2),[],1),[],2)-...
                min(min(point.profile(:,:,2),[],1),[],2);
            om1=point.omega1;
            om2=point.omega2;
            r1=(deltax*om2-deltaxdot)/2/(om2-om1);
            r2=(deltax*om1-deltaxdot)/2/(om1-om2);
            r=r1*sqrt(2*(1-cos(om1*tau)))+r2*sqrt(2*(1-cos(om2*tau)));
            pbl(k,2)=r;
        end
    end
    handle=plot(pbl(:,1),pbl(:,2),clr);
    
    
    
    
end
return
    function ukl=j2kl(varargin)
        if length(varargin)==1
            ukl=permute(reshape(varargin{1},[n, Nom1, Nom2]),[2 3 1]);
        elseif length(varargin)==2
            ukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{2}]),[2 3 1]);
        elseif length(varargin)==3
            ukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{3}]),[2 3 1]);
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
            n1=varargin{1};
            n2=varargin{2};
        end
        ukl=varargin{1};
        uj=reshape(permute(ukl,[3 1 2]),[n1*n2*size(ukl,3) 1]);
    end
end

