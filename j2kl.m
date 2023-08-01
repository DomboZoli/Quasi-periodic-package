    function pukl=j2kl(varargin)
%%  COPYRIGHT
% Zoltan Dombovari, dombovari@mm.bme.hu, 
% Department of Applied Mechanics, 
% Faculty of Mechanical Engineering
% Budapest University of Technology and Economics
% statement: This is a purely research oriented algortihm, made in a result oriented manner. It is only optimized up to a convenient level. I apologise all inefficiency, errors and grammatic mistakes and lack in/of comments. Please report suggestions on the above email. Any use or publications based on the algorithm must be authorized by the author 
% optimized for matlab 2018b
%% main code: outsources due to parallellization
    n=4;
        if length(varargin)==1
            pukl=permute(reshape(varargin{1},[n, Nom1, Nom2]),[2 3 1]);
        elseif length(varargin)==2
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{2}]),[2 3 1]);
        elseif length(varargin)==3
            pukl=permute(reshape(varargin{1},[n, varargin{2}, varargin{3}]),[2 3 1]);
        elseif length(varargin)==4
            pukl=permute(reshape(varargin{1},[varargin{4}, varargin{2}, varargin{3}]),[2 3 1]);
        end
    end