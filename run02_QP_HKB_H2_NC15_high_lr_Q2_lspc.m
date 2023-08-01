%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% SCRIPT
load('br_ahigh_left_QP2_NC31_lspc');

%% left
"br_ahigh_left_QP2_NC15_lspc"
br_ahigh_left_QP2_NC15_lspc=br_ahigh_left_QP2_NC31_lspc;
br_ahigh_left_QP2_NC15_lspc.points=br_ahigh_left_QP2_NC31_lspc.points(1);
br_ahigh_left_QP2_NC15_lspc.points.num.Nc=round(1.2*numel(br_ahigh_left_QP2_NC31_lspc.points));
% br_ahigh_left_QP1_NC31_lspc.points.par.num.ds=-br_ahigh_left_QP1_NC31_lspc.points.par.num.ds;
par=br_ahigh_left_QP2_NC15_lspc.points(1).par;
par.num.Nc=25;
par.num.pctype=1; %LS phasecondition
par.num.scheme=2; %Legendre
par.num.NC1=15;
par.num.NC2=15;
%originally 31
par.num.CONV=0;
par.num.ds=5*br_ahigh_left_QP2_NC31_lspc.points(1).par.num.ds;
par.num.stepsavename='bkpbranch_QP_HKB_H2_NC15_high_left_Q2_lspc';
try
    br_ahigh_left_QP2_NC15_lspc=cont_2dtori(br_ahigh_left_QP2_NC15_lspc.points(end),par,br_ahigh_left_QP2_NC15_lspc.contpar);
%     br_ahigh_left_QP2_NC31_lspc.points(end).par.num.Nc=50;br_ahigh_left_QP2_NC31_lspc=cont_2dtori(br_ahigh_left_QP2_NC31_lspc);
catch
    global bkpbranch;
    br_ahigh_left_QP2_NC15_lspc=bkpbranch;
end
save('br_ahigh_left_QP2_NC15_lspc.mat','br_ahigh_left_QP2_NC15_lspc');    
    

%% right %this breaks
"br_ahigh_right_QP2_NC15_lspc"
br_ahigh_right_QP2_NC15_lspc=br_ahigh_left_QP2_NC31_lspc;
br_ahigh_right_QP2_NC15_lspc.points=br_ahigh_left_QP2_NC31_lspc.points(1); 
br_ahigh_right_QP2_NC15_lspc.points.num.Nc=round(1.2*numel(br_ahigh_left_QP2_NC31_lspc.points));
% br_ahigh_left_QP1_NC31_lspc.points.par.num.ds=-br_ahigh_left_QP1_NC31_lspc.points.par.num.ds;
par=br_ahigh_right_QP2_NC15_lspc.points(1).par;
par.num.Nc=50;
par.num.pctype=1; %LS phasecondition
par.num.scheme=2; %Legendre
par.num.NC1=15;
par.num.NC2=15;
%originally 31
par.num.CONV=0;
par.num.ds=5*br_ahigh_left_QP2_NC31_lspc.points(1).par.num.ds;
par.num.stepsavename='bkpbranch_QP_HKB_H2_NC15_high_right_Q2_lspc';
try
    br_ahigh_right_QP2_NC15_lspc=cont_2dtori(br_ahigh_right_QP2_NC15_lspc.points(end),par,br_ahigh_right_QP2_NC15_lspc.contpar);
%     br_ahigh_left_QP2_NC31_lspc.points(end).par.num.Nc=50;br_ahigh_left_QP2_NC31_lspc=cont_2dtori(br_ahigh_left_QP2_NC31_lspc);
catch
    global bkpbranch;
    br_ahigh_right_QP2_NC15_lspc=bkpbranch;
end
save('br_ahigh_right_QP2_NC15_lspc.mat','br_ahigh_right_QP2_NC15_lspc');  
    