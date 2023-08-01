function f=sys_rhs(xx,par)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
%% main code: original DDE-BIFTOOL file

% alpha beta gamma omega a b tau
%         1:x_1,2: x_2, 3:dx_1, 4:dx_2   


f(1,1)=xx(3,1);
f(2,1)=xx(4,1);
f(3,1)=-(xx(3,1).*(par(1).*xx(1,1).^2+par(2).*xx(3,1).^2-par(3))+par(4).^2.*xx(1,1))...
    +(par(5)+par(6).*(xx(1,1)-xx(2,2)).^2).*(xx(3,1)-xx(4,2));
f(4,1)=-(xx(4,1).*(par(1).*xx(2,1).^2+par(2).*xx(4,1).^2-par(3))+par(4).^2.*xx(2,1))...
    + (par(5)+par(6).*(xx(2,1)-xx(1,2)).^2).*(xx(4,1)-xx(3,2));

return;
