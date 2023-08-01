%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% SCRIPT
%script switch to higher resolution on QP_HKB_H2_NC25_Q1
load('br_a_left_QP1');
%% left
"br_a_left_QP1_lspc"
br_a_left_QP1_lspc=br_a_left_QP1;
br_a_left_QP1_lspc.points=br_a_left_QP1.points(1);
br_a_left_QP1_lspc.points.Nc=round(1.2*numel(br_a_left_QP1.points));
par=br_a_left_QP1_lspc.points(1).par;
par.num.Nc=30;
par.num.pctype=1; %LS phasecondition
par.num.scheme=2; %Legendre
par.num.NC1=15;
par.num.NC2=15;
par.num.CONV=0;
par.num.ds=-5*br_a_left_QP1_lspc.points.par.num.ds;
par.num.stepsavename='bkpbranch_br_a_left_QP1_NC31_lspc';
try
    br_a_left_QP1_lspc=cont_2dtori(br_a_left_QP1_lspc.points(end),par,br_a_left_QP1_lspc.contpar);
catch
    global bkpbranch;
    br_a_left_QP1_lspc=bkpbranch;
end
save('br_a_left_QP1_lspc.mat','br_a_left_QP1_lspc');    
%% right
 "br_a_right_QP1_lspc"
br_a_right_QP1_lspc=br_a_left_QP1;
br_a_right_QP1_lspc.points=br_a_left_QP1.points(1);
br_a_right_QP1_lspc.points.Nc=round(1.2*numel(br_a_left_QP1.points));
par=br_a_right_QP1_lspc.points(1).par;
par.num.Nc=35;
par.num.pctype=1; %LS phasecondition
par.num.scheme=2; %Legendre
par.num.NC1=15;
par.num.NC2=15;
par.num.CONV=0;
par.num.ds=5*br_a_right_QP1_lspc.points.par.num.ds;
par.num.stepsavename='bkpbranch_br_a_left_QP1_NC31_lspc';
try
    br_a_right_QP1_lspc=cont_2dtori(br_a_right_QP1_lspc.points(end),par,br_a_right_QP1_lspc.contpar);
catch
    global bkpbranch;
    br_a_right_QP1_lspc=bkpbranch;
end
save('br_a_right_QP1_lspc.mat','br_a_right_QP1_lspc');
