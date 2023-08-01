function vcijtau=cptaustate_ex(vkl,s1,s2,inds1,inds2,numpars,deriv)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code:outsourced for parallellization
%not nested cptaustate function for parfor
%numpar numparameters
Nom1=numpars.Nom1;
Nom2=numpars.Nom2;
P=numpars.P;
n=numpars.n;
%length(inds1)=P+1
%length(inds2)=P+1
%part delayed approx
%extension part
vkld=cat(2,vkl(1:Nom1-1,1:Nom2-1,:),vkl(1:Nom1-1,1:Nom2,:));
vkldd=cat(1,vkld,vkld(1:Nom1-1,:,:));
vkldd(2*Nom1-1,:,:)=vkldd(1,:,:);
vkldd(:,2*Nom2-1,:)=vkldd(:,1,:);
%         vkldd
%the delay is produced in the representation points and interpolated at the original collocation points
phis1=linspace(-2*pi,2*pi,Nom1*2-1);
phis2=linspace(-2*pi,2*pi,Nom2*2-1);
vkltau=zeros(P+1,P+1,n);
for l_int=1:n
    vkltau(:,:,l_int)=interp2(phis1,phis2,vkldd(:,:,l_int)',phis1(Nom1-1+inds1)-mod(s1,2*pi),(phis2(Nom2-1+inds2)-mod(s2,2*pi)).','linear')';
end
%part states at the collocation points
vcijtau=cpstate_ex(vkltau,numpars,deriv);

