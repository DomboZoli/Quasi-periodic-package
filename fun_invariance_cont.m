function ices=fun_invariance_cont(pj,par,cpari,upj,contstepcond,numpars)
% continuation, condition and the invariance equation combined for quasiperiod vibration
%% INPUT
%   pj: actual (jth) iteration step, [puij; pom;]
%       puij: (p: not iterated) representation points 
%       pom: (p: not iterated) vibration frequency
%   par: parameterlist (see parameterlist file)
%   contstepcond: 0: pseudo arclengh, 1: constant np parameter
%% OUTPUT
%   jc: colocation Jacobian
%% VERSION v2.5
%   v1.0: not optimized brute force solution to have a solution
%   v1.0b: modified for 1dof turning
%   v2.0: fixed frequency truncation by par.num.fixom
%   v2.1: distance debugged
%   v2.2: intruducing constant parameter continuation condition
%   v2.3: dim and dof in the par.system list, def_taus introduced
%   v2.4: numpars has introduced for external jac_cond
%   v2.5: finalized
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% preset parameters
if ~isempty(numpars)
    %collocation order
    P=par.num.p;
    %collocation intervals for 2D tori both directions
    NC1=par.num.NC1;
    NC2=par.num.NC2;
    n=2*par.system.dof;
    %number of parameters
    np=length(cpari);
else
    P=numpars.P;
    n=numpars.n;
    NC1=numpars.NC1;
    NC2=numpars.NC2;
    np=numpars.np;
end
%number 
%% main code
ices=zeros((NC1*P+1)*(NC2*P+1)*n+2+np,1);
% include invariance equation
ices(1:(NC1*P+1)*(NC2*P+1)*n+2,1)=fun_invariance_2dtori(pj,par);
%parameter continuations and conditions
if (np==1 && isempty(upj)) || (np==2)
    %store the minimizing value
    ices((NC1*P+1)*(NC2*P+1)*n+2+1)=par.system.sys_cond(pj,par,[],numpars);
elseif ~isempty(upj)
    %distance condition for continuation
    %     ices((NC1*P+1)*(NC2*P+1)*n+2+np)=upj.'*upj-2*upj.'*pj+pj.'*pj-par.num.ds^2;
    switch contstepcond
        case 0
            %pseudo arc length
            ices((NC1*P+1)*(NC2*P+1)*n+2+np)=(upj-pj).'*(upj-pj)-par.num.ds^2;
        case 1
            %constant parameter
            ices((NC1*P+1)*(NC2*P+1)*n+2+np)=pj((NC1*P+1)*(NC2*P+1)*n+2+np)-upj((NC1*P+1)*(NC2*P+1)*n+2+np);
    end
end
%truncation
if par.num.fixom>0
    ices=[ices(1:(NC1*P+1)*(NC2*P+1)*n+par.num.fixom-1,1);ices((NC1*P+1)*(NC2*P+1)*n+par.num.fixom+1:end,1);];
end

return