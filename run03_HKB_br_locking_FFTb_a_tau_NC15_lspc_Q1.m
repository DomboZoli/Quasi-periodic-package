%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% SCRIPT
%name
name='br_a_tau_left_QP1_NC15_lspc';
%plotbranches
load('br_a_right_QP1_lspc.mat');
load('numpars_NC31.mat');
% % 
%% left
%initial point that has sign of break
STP=25;
disp('bkpbranch_br_a_tau_left_NC15_Q2_lspc');
par=br_a_right_QP1_lspc.points(STP).par;
par.num.maxNRIc=20;
par.num.absTOL=1e-5;
par.num.parTOL=1e-5*ones(size(par.num.parTOL));
par.num.ds=1.5/5;
par.num.Nc=100;
par.num.NC1=15;
par.num.NC2=15;
%the higher harmonics relative error on the torus: crucial to check by test run breaked in fun_rhs_HKB to have epsilon value
par.system.Dom=sqrt(0.00528466564013678);% using fun_cond_HKB_FFT to extract epsilon
%set harmonics range
par.system.condpar.iTTwmin=2;
par.system.condpar.iTTwmax=15;
% par.system.Dom=0.8706
% par.system.Dom=res_a_left.points(122).omega1-res_a_left.points(122).omega2;
par.parlist={'alp'  'beta'  'gam'  'om'  'a'  'b'  'tau' 'Dom'};
par.num.pctype=1;
br_a_tau_left_NC15_Q1_lspc=br_a_right_QP1_lspc;
br_a_tau_left_NC15_Q1_lspc.contpar=[5 7];
br_a_tau_left_NC15_Q1_lspc.points=br_a_right_QP1_lspc.points(STP);
br_a_tau_left_NC15_Q1_lspc.points.par=par;
br_a_tau_left_NC15_Q1_lspc.points.par.system.sys_cond=@fun_cond_HKB_FFT;
par=br_a_tau_left_NC15_Q1_lspc.points(1).par;
par.num.stepsavename='bkpbranch_br_a_tau_left_NC15_Q1_lspc';
% par.system.Dom=0.981;
try
    br_a_tau_left_NC15_Q1_lspc=cont_2dtori(br_a_tau_left_NC15_Q1_lspc.points(1),par,br_a_tau_left_NC15_Q1_lspc.contpar);
% br_a_tau_left_NC15_Q1_lspc=cont_2dtori(br_a_tau_left_NC15_Q1_lspc);
catch
    global bkpbranch;
    br_a_tau_left_NC15_Q1_lspc=bkpbranch;
end
save('br_a_tau_left_NC15_Q1_lspc.mat','br_a_tau_left_NC15_Q1_lspc');

%% right
STP=25;
disp('bkpbranch_br_a_tau_left_NC15_Q2_lspc');
par=br_a_right_QP1_lspc.points(STP).par;
par.num.maxNRIc=20;
par.num.absTOL=1e-5;
par.num.parTOL=1e-5*ones(size(par.num.parTOL));
par.num.ds=1.5/5;
par.num.Nc=100;
par.num.NC1=15;
par.num.NC2=15;
%the higher harmonics relative error on the torus
par.system.Dom=sqrt(0.00528466564013678);%
%set harmonics range
par.system.condpar.iTTwmin=2;
par.system.condpar.iTTwmax=15;
% par.system.Dom=0.8706
% par.system.Dom=res_a_left.points(122).omega1-res_a_left.points(122).omega2;
par.parlist={'alp'  'beta'  'gam'  'om'  'a'  'b'  'tau' 'Dom'};
par.num.pctype=1;
br_a_tau_right_NC15_Q1_lspc=br_a_right_QP1_lspc;
br_a_tau_right_NC15_Q1_lspc.contpar=[5 7];
br_a_tau_right_NC15_Q1_lspc.points=br_a_right_QP1_lspc.points(STP);
br_a_tau_right_NC15_Q1_lspc.points.par=par;
br_a_tau_right_NC15_Q1_lspc.points.par.system.sys_cond=@fun_cond_HKB_FFT;
par=br_a_tau_right_NC15_Q1_lspc.points(1).par;
par.num.stepsavename='bkpbranch_br_a_tau_right_NC15_Q1_lspc';
% par.system.Dom=0.981;
try
    br_a_tau_right_NC15_Q1_lspc=cont_2dtori(br_a_tau_right_NC15_Q1_lspc.points(1),par,br_a_tau_right_NC15_Q1_lspc.contpar);
catch
    global bkpbranch;
    br_a_tau_right_NC15_Q1_lspc=bkpbranch;
end
save('br_a_tau_right_NC15_Q1_lspc.mat','br_a_tau_right_NC15_Q1_lspc');
