function res=fun_deri_HKB(ph1,ph2,y,ytau,par,ns)
%derivatives for HKB definition, first order vectorfield
%% INPUT
%   ph1: optional phase 1, 
%   ph2: optional phase 2, 
%   y: state size(y)=[n 1],
%   ytau: delayes states size(ytau)=[n Ntau],
%   par: parameter list (see parameterlist file)
%   ns: switch
%% OUTPUT
%   res: colocation Jacobian
%% VERSION v1.0
%   v1.0: not optimized brute force solution to have a solution
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code
%FM
FM=1;
%system DOF
dof=2;
% alpha beta gamma omega a b tau
parbt=[par.system.alp par.system.beta par.system.gam par.system.om par.system.a par.system.b par.system.tau];

if ns==1
    %exact derivative w.r.t. the present state
    res=sys_deri([y ytau],parbt,0,[],[]);
elseif ns>1 && ns<=FM+1
    %derivative w.r.t. the delayed state
    res=sys_deri([y ytau],parbt,1,[],[]);
elseif ns==FM+1+1
    %derivative w.r.t. the first parameter "alpha"
    res=sys_deri([y ytau],parbt,[],1,[]);
elseif ns==FM+1+2
    %derivative w.r.t. "beta"
    res=sys_deri([y ytau],parbt,[],2,[]);
elseif ns==FM+1+3
    %derivative w.r.t. "gamma"
    res=sys_deri([y ytau],parbt,[],3,[]);
elseif ns==FM+1+4
    %derivative w.r.t. "omega"
    res=sys_deri([y ytau],parbt,[],4,[]);
elseif ns==FM+1+5
    %derivative w.r.t. "a"
    res=sys_deri([y ytau],parbt,[],5,[]);
elseif ns==FM+1+6
    %derivative w.r.t. "b"
    res=sys_deri([y ytau],parbt,[],6,[]);
elseif ns==FM+1+7
    %derivative w.r.t. "tau"
    res=sys_deri([y ytau],parbt,[],7,[]);
elseif ns==1i
    %partial om1 related explicit time
    res=zeros(2*dof,1);
elseif ns==2i
    %partial om2 related explicit time
    res=zeros(2*dof,1);
end

return

end