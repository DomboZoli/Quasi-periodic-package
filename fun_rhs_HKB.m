function fcij=fun_rhs_HKB(ucij,ucijtau,c1ij,par,om1,om2)
% HKB definition, first order vectorfield
%% INPUT
%   ucij: actual (jth) state, 
%   ucijtau: delayed state
%   par: parameter list (see parameterlist file)
%   om1: freq1
%   om2: freq2
%   TODO: ucijtau multidelay!!!
%% OUTPUT
%   jc: colocation Jacobian
%% VERSION v1.1
%   v1.0: not optimized brute force solution to have a solution
%   v1.1: preallocation fcij
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% preset parameters
%including flyover map

%define parameters for sys_rhs
% alpha beta gamma omega a b tau
parbt=[par.system.alp par.system.beta par.system.gam par.system.om par.system.a par.system.b par.system.tau];
%preallocate fcij
fcij=zeros(size(c1ij,1),size(c1ij,2),size(ucij,3));
%brute force method
for i=1:size(c1ij,1)
    for j=1:size(c1ij,2)
        %rhs definition, not really an efficient one
%         1:x_1,2: x_2, 3:dx_1, 4:dx_2    
        fcij(i,j,:)=sys_rhs([permute(ucij(i,j,:),[3 2 1]) permute(ucijtau(i,j,:),[3 2 1])],parbt);
    end
end

return
end
