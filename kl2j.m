    function uj=kl2j(varargin)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code: outsourced due to parallellization
            pukl=varargin{1};
        if length(varargin)==1
            n1=size(pukl,1);
            n2=size(pukl,2);
        elseif length(varargin)==2
            n1=varargin{2};
            n2=n1;
        else
            n1=varargin{2};
            n2=varargin{3};
        end;

        uj=reshape(permute(pukl,[3 1 2]),[n1*n2*size(pukl,3) 1]);
    end
