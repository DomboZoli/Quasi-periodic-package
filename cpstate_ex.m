function pvcij=cpstate_ex(pvkl,numpars,deriv)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code
%numpar numparameters
P=numpars.P;
n=numpars.n;
cPpq=numpars.cPpq;
cd1Ppq=numpars.cd1Ppq;
cd2Ppq=numpars.cd2Ppq;
%part approximation, size(pvkl)=[P+1 P+1 n], size(pcvkl)=[P P n]
%declarate the result function
%states at the collocation points
%part of the representation points
switch deriv
    case 0
        pvcij=permute(sum(sum(repmat(permute(pvkl(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
            .*repmat(cPpq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
    case 1
        pvcij=permute(sum(sum(repmat(permute(pvkl(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
            .*repmat(cd1Ppq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
        
    case 2
        pvcij=permute(sum(sum(repmat(permute(pvkl(1:P+1,1:P+1,:),[4 5 1 2 3]),[P,P,1,1,1])...
            .*repmat(cd2Ppq,[1 1 1 1 n]),3),4),[1 2 5 3 4]);
end
