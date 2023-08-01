%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% scripts
%plotbranches
load('br_ahigh_left_QP2_NC15_lspc.mat');
load('numpars_NC31.mat');

%% left
STP=20;
par=br_ahigh_left_QP2_NC15_lspc.points(STP).par;

% par.system.Dom=res_a_left.points(122).omega1-res_a_left.points(122).omega2;
par.parlist={'alp'  'beta'  'gam'  'om'  'a'  'b'  'tau' 'Dom'};

br_ahigh_tau_left_NC15_lspc=br_ahigh_left_QP2_NC15_lspc;
br_ahigh_tau_left_NC15_lspc.contpar=[5 7];
br_ahigh_tau_left_NC15_lspc.points=br_ahigh_left_QP2_NC15_lspc.points(STP);
br_ahigh_tau_left_NC15_lspc.points.par=par;
br_ahigh_tau_left_NC15_lspc.points.par.system.sys_cond=@fun_cond_HKB_FFT;
par=br_ahigh_tau_left_NC15_lspc.points(1).par;
par.num.maxNRIc=40;
par.num.maxNRc=80;
par.num.absTOL=1e-5;
par.num.parTOL=1e-5*ones(size(par.num.parTOL));
par.num.ds=1.5/3;
par.num.Nc=50;
par.num.NC1=15;
par.num.NC2=15;
%the higher harmonics relative error on the torus: crucial to check by test run breaked in fun_rhs_HKB to have epsilon value
par.system.Dom=sqrt(0.0135882767135607); %%% using fun_cond_HKB_FFT to extract epsilon
%set harmonics range
par.system.condpar.iTTwmin=2;
par.system.condpar.iTTwmax=15;
par.num.stepsavename='bkpbranch_br_ahigh_tau_left_NC15_lspc';
try
    br_ahigh_tau_left_NC15_lspc=cont_2dtori(br_ahigh_tau_left_NC15_lspc.points(1),par,br_ahigh_tau_left_NC15_lspc.contpar);
catch
    global bkpbranch;
    br_ahigh_tau_left_NC31_lspc_P8_26=bkpbranch;
end
save('br_ahigh_tau_left_NC15_lspc.mat','br_ahigh_tau_left_NC15_lspc');

%% right

STP=20;
par=br_ahigh_left_QP2_NC15_lspc.points(STP).par;

% par.system.Dom=res_a_left.points(122).omega1-res_a_left.points(122).omega2;
par.parlist={'alp'  'beta'  'gam'  'om'  'a'  'b'  'tau' 'Dom'};

br_ahigh_tau_right_NC15_lspc=br_ahigh_left_QP2_NC15_lspc;
br_ahigh_tau_right_NC15_lspc.contpar=[5 7];
br_ahigh_tau_right_NC15_lspc.points=br_ahigh_left_QP2_NC15_lspc.points(STP);
br_ahigh_tau_right_NC15_lspc.points.par=par;
br_ahigh_tau_right_NC15_lspc.points.par.system.sys_cond=@fun_cond_HKB_FFT;
par=br_ahigh_tau_right_NC15_lspc.points(1).par;
par.num.maxNRIc=40;
par.num.maxNRc=80;
par.num.absTOL=1e-5;
par.num.parTOL=1e-5*ones(size(par.num.parTOL));
par.num.ds=-1.5/3;
par.num.Nc=50;
par.num.NC1=15;
par.num.NC2=15;
%the higher harmonics relative error on the torus: crucial to check by test run breaked in fun_rhs_HKB to have epsilon value
par.system.Dom=sqrt(0.0135882767135607); %%fun_cond_HKB_FFT v11.2
%set harmonics range
par.system.condpar.iTTwmin=2;
par.system.condpar.iTTwmax=15;
par.num.stepsavename='bkpbranch_br_ahigh_tau_right_NC15_lspc';
try
    br_ahigh_tau_right_NC15_lspc=cont_2dtori(br_ahigh_tau_right_NC15_lspc.points(1),par,br_ahigh_tau_right_NC15_lspc.contpar);
catch
    global bkpbranch;
    br_ahigh_tau_left_NC31_lspc_P8_26=bkpbranch;
end
save('br_ahigh_tau_right_NC15_lspc.mat','br_ahigh_tau_right_NC15_lspc');
