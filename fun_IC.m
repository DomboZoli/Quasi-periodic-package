%linear indexes for collocation vector (i,j) interval
function res=fun_IC(i,j,numpars)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code: outsourced due to parallel algorithm
P=numpars.P;
n=numpars.n;
NC1=floor((numpars.Nom1-1)/P);
%         NC2=floor((numpars.Nom2-1)/P);
res=reshape(repmat(((j-1)*P+(1:P)-1)*NC1*P*n,[P 1 n])+repmat(permute(((i-1)*P+(1:P)-1)*n,[2 1 3]),[1 P n])+repmat(permute(1:n,[1 3 2]),[P P 1]),[P^2*n 1]);
end