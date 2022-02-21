function res = crop(varargin)
% Function crop is to crop 1d, 2d  and 3d data around the center
%
% Input Variables:
% varagin - original data and central region size
%
% Output Variables:
% res - data after croped
%
% Example:
% 1D vector x: res = crop(x,sx);
% 2D matrix x: res = crop(x,sx,sy);
% 3D matrix x: res = crop(x,sx,sy,sz);
%
% Record of Revisions:
% Jun-13-2020===Zhao He===add 1D and 3D data processing functions based on
%                           'crop2D.m' writen by Michael Lustig 2007.
% Jun-17-2020===Yuan Feng===add comments
% Jun-17-2020===Zhao He===allow the odd size of input data and central roi

x = varargin{1};

if isvector(x)
    if (diff(size(x))>0) % if a row vector
        mx = size(x,2);
        sx = varargin{2};
        idxx = floor(mx/2)+1 - floor(sx/2) : floor(mx/2) + ceil(sx/2);
        res = x(1,idxx);
    else % if a column vector
        mx = size(x,1);
        sx = varargin{2};
        idxx = floor(mx/2)+1 - floor(sx/2) : floor(mx/2) + ceil(sx/2);
        res = x(idxx,1);
    end
    
else
    
    n = ndims(x); 
    if n==2 % if 2d data
        [mx,my] = size(x);
        sx = varargin{2}; sy = varargin{3};
        idxx = floor(mx/2)+1 - floor(sx/2) : floor(mx/2) + ceil(sx/2);
        idyy = floor(my/2)+1 - floor(sy/2) : floor(my/2) + ceil(sy/2);
        res = x(idxx,idyy);
    elseif n==3 % if 3d data
        [mx,my,mz] = size(x);
        sx = varargin{2}; sy = varargin{3}; sz = varargin{4};
        idxx = floor(mx/2)+1 - floor(sx/2) : floor(mx/2) + ceil(sx/2);
        idyy = floor(my/2)+1 - floor(sy/2) : floor(my/2) + ceil(sy/2);
        idzz = floor(mz/2)+1 - floor(sz/2) : floor(mz/2) + ceil(sz/2);
        res = x(idxx,idyy,idzz);
    else
        disp('Input must be 1D, 2D or 3D data!');return;
    end
end






