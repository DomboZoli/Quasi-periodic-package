%linear indexes for representation vector at (k,l) interval
function res=fun_IR(k,l,numpars)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code: outsourced due to parallellization 
P=numpars.P;
n=numpars.n;
NC1=floor((numpars.Nom1-1)/P);
res=reshape(repmat(((l-1)*P+(1:P+1)-1)*(NC1*P+1)*n,[P+1 1 n])+repmat(permute(((k-1)*P+(1:P+1)-1)*n,[2 1 3]),[1 P+1 n])+repmat(permute(1:n,[1 3 2]),[P+1 P+1 1]),[(P+1)^2*n 1]);
end