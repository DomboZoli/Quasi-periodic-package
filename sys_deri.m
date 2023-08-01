function J=sys_deri(xx,par,nx,np,v)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
%% original DDE-BIFTOOL file
% alpha beta gamma omega a b tau

J=[];

if length(nx)==1 & length(np)==0 & isempty(v)
    % first order derivatives wrt state variables
    if nx==0 % derivative wrt x(t)
        J(1,3)=1;
        J(2,4)=1;
        J(3,1)=-2.*par(1).*xx(3,1).*xx(1,1)-par(4).^2+...
            2.*par(6).*(xx(1,1)-xx(2,2)).*(xx(3,1)-xx(4,2));
        J(3,3)=-par(1).*xx(1,1).^2-3.*par(2).*xx(3,1).^2+par(3)+...
            par(5)+par(6).*(xx(1,1)-xx(2,2)).^2;
        J(4,2)=-2.*par(1).*xx(4,1).*xx(2,1)-par(4).^2+...
            2.*par(6).*(xx(2,1)-xx(1,2)).*(xx(4,1)-xx(3,2));
        J(4,4)=-par(1).*xx(2,1).^2-3.*par(2).*xx(4,1).^2+par(3)+...
            par(5)+par(6).*(xx(2,1)-xx(1,2)).^2;
    elseif nx==1 %derivative wrt x(t-tau)
        J(3,2)=-2.*par(6).*(xx(1,1)-xx(2,2)).*(xx(3,1)-xx(4,2));
        J(3,4)=-(par(5)+par(6).*(xx(1,1)-xx(2,2)).^2);
        J(4,1)=-2.*par(6).*(xx(2,1)-xx(1,2)).*(xx(4,1)-xx(3,2));
        J(4,3)=-(par(5)+par(6).*(xx(2,1)-xx(1,2)).^2);
    end
elseif length(nx)==0 & length(np)==1 & isempty(v)
    %first order derivatives wrt parameters
    if np==1 %derivative wrt alpha
        J(3,1)=-xx(3,1).*xx(1,1).^2;
        J(4,1)=-xx(4,1).*xx(2,1).^2;
    elseif np==2 %derivative wrt beta
        J(3,1)=-xx(3,1).^3;
        J(4,1)=-xx(4,1).^3;
    elseif np==3 %derivative wrt gamma
        J(3,1)=xx(3,1);
        J(4,1)=xx(4,1);
    elseif np==4 %derivative wrt omega
        J(3,1)=-2.*par(4).*xx(1,1);
        J(4,1)=-2.*par(4).*xx(2,1);
    elseif np==5 %derivative wrt a
        J(3,1)=xx(3,1)-xx(4,2);
        J(4,1)=xx(4,1)-xx(3,2);
    elseif np==6 %derivative wrt b
        J(3,1)=(xx(1,1)-xx(2,2)).^2.*(xx(3,1)-xx(4,2));
        J(4,1)=(xx(2,1)-xx(1,2)).^2.*(xx(4,1)-xx(3,2));
    elseif np==7 %derivative wrt tau
        J=zeros(4,1);
    end
elseif length(nx)==1 & length(np)==1 & isempty(v)
    %mixed state, parameter derivatives
    if nx==0 %derivative wrt x(t)
        if np==1 %derivative wrt alpha
            J(3,1)=-2.*xx(3,1).*xx(1,1);
            J(3,3)=-xx(1,1).^2;
            J(4,2)=-2.*xx(4,1).*xx(2,1);
            J(4,4)=-xx(2,1).^2;
        elseif np==2
            J(3,3)=-3.*xx(3,1).^2;
            J(4,4)=-3.*xx(4,1).^2;
        elseif np==3
            J(3,3)=1;
            J(4,4)=1;
        elseif np==4
            J(3,1)=-2.*par(4);
            J(4,2)=-2.*par(4);
            J(4,4)=0;
        elseif np==5
            J(3,3)=1;
            J(4,4)=1;
        elseif np==6
            J(3,1)=2.*(xx(1,1)-xx(2,2)).*(xx(3,1)-xx(4,2));
            J(3,3)=(xx(1,1)-xx(2,2)).^2;
            J(4,2)=2.*(xx(2,1)-xx(1,2)).*(xx(4,1)-xx(3,2));
            J(4,4)=(xx(2,1)-xx(1,2)).^2;
        elseif np==7
            J=zeros(4,4);
        end
    elseif nx==1
        if np==5
            J(3,4)=-1;
            J(4,3)=-1;
        elseif np==6
            J(3,2)=-2.*(xx(1,1)-xx(2,2)).*(xx(3,1)-xx(4,2));
            J(3,4)=-(xx(1,1)-xx(2,2)).^2;
            J(4,1)=-2.*(xx(2,1)-xx(1,2)).*(xx(4,1)-xx(3,2));
            J(4,3)=-(xx(2,1)-xx(1,2)).^2;
        else
            J=zeros(4,4);
        end
    end
elseif length(nx)==2 & length(np)==0 & ~isempty(v)
    %second order derivatives wrt state variables
    if nx(1)==0
        if nx(2)==0
            J(3,1)=v(1).*(-2.*par(1).*xx(3,1)+2.*par(6).*(xx(3,1)-xx(4,2)))+...
                v(3).*(-2.*par(1).*xx(1,1)+2.*par(6).*(xx(1,1)-xx(2,2)));
            J(3,3)=v(1).*(-2.*par(1).*xx(1,1)+2.*par(6).*(xx(1,1)-xx(2,2)))+...
                v(3).*(-6.*par(2).*xx(3,1));
            J(4,2)=v(2).*(-2.*par(1).*xx(4,1)+2.*par(6).*(xx(4,1)-xx(3,2)))+...
                v(4).*(-2.*par(1).*xx(2,1)+2.*par(6).*(xx(2,1)-xx(1,2)));
            J(4,4)=v(2).*(-2.*par(1).*xx(2,1)+2.*par(6).*(xx(2,1)-xx(1,2)))+...
                v(4).*(-6.*par(2).*xx(4,1));
        elseif nx(2)==1
            J(3,2)=v(1).*(-2.*par(6).*(xx(3,1)-xx(4,2)))-v(3).*2.*par(6).*(xx(1,1)-xx(2,2));
            J(3,4)=v(1).*(-2.*par(6).*(xx(1,1)-xx(2,2)));
            J(4,1)=v(2).*(-2.*par(6).*(xx(4,1)-xx(3,2)))-v(4).*2.*par(6).*(xx(2,1)-xx(1,2));
            J(4,3)=v(2).*(-2.*par(6).*(xx(2,1)-xx(1,2)));
        end
    elseif nx(1)==1
        if nx(2)==0
            J(3,1)=v(2).*(-2.*par(6).*(xx(3,1)-xx(4,2)))-v(4).*2.*par(6).*(xx(1,1)-xx(2,2));
            J(3,3)=v(2).*(-2.*par(6).*(xx(1,1)-xx(2,2)));
            J(4,2)=v(1).*(-2.*par(6).*(xx(4,1)-xx(3,2)))-v(3).*2.*par(6).*(xx(2,1)-xx(1,2));
            J(4,4)=v(1).*(-2.*par(6).*(xx(2,1)-xx(1,2)));
        elseif nx(2)==1
            J(3,2)=-v(2).*(-2.*par(6).*(xx(3,1)-xx(4,2)))-v(4).*2.*par(6).*(xx(1,1)-xx(2,2));
            J(3,4)=-v(2).*(-2.*par(6).*(xx(1,1)-xx(2,2)));
            J(4,1)=-v(1).*(-2.*par(6).*(xx(4,1)-xx(3,2)))-v(3).*2.*par(6).*(xx(2,1)-xx(1,2));
            J(4,3)=-v(1).*(-2.*par(6).*(xx(2,1)-xx(1,2)));
        end
    end
            

end
        
if isempty(J)
  err=[nx np size(v)]
  error('SYS_DERI: requested derivative could not be computed!');
end;


return;