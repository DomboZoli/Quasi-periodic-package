function [taus,npt]=fun_sys_tau_HKB(pj,par)
%definition of delayes for the milling model
% work after the version >v9.2 for cont_2dtori
%npt: delay parameter index
%% VERSIONS 1.1
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code
%WARNING NFO==1  delay
taus=par.system.tau;
npt=7;
end