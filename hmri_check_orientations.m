function [sts, str] = hmri_check_orientations(V, verbose)
% Check the dimensions and orientations of the images
% FORMAT [sts, str] = hmri_check_orientations(V [,verbose])
% V       - a struct array as returned by spm_vol
% verbose - [Default: true]
%
% sts     - status (true means OK)
% str     - string describing status, empty if OK
%
% When used without LHS, this function throws an error accordingly.
%__________________________________________________________________________
% Adapted from John Ashburner's spm_check_orientations.m
% to save all console output in str
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

sts = true;
str = '';

if nargin < 2, verbose = true; end

dims = cat(1,V.dim);
if any(any(diff(dims,1,1),1))
    sts = false;
    str = sprintf('%s\n\t** The images do not all have the same dimensions **\n',str);
    if verbose
        str = sprintf('%sThe function assumes that a voxel in one image corresponds with\n',str);
        str = sprintf('%sthe same  voxel in another.   This  is not a safe assumption if\n',str);
        str = sprintf('%sthe  image dimensions differ.   Please  ensure  that  you  have\n',str);
        str = sprintf('%sprocessed all the image data in the same way (eg. check spatial\n',str);
        str = sprintf('%snormalisation bounding-boxes, voxel-sizes etc).\n',str);
        str = sprintf('%sHere are the dimensions of the image volumes.  This list can be\n',str);
        str = sprintf('%sused to determine which file(s) are causing the problem.\n\n',str);
        for i=1:numel(V)
            str = sprintf('%s[%d %d %d]  %s\n', str, V(i).dim, V(i).fname);
        end
        str = sprintf('%s\n',str);
    end
    if ~nargout, error('The dimensions must be identical for this procedure.'); end
end

matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4))
    sts = false;
    str = sprintf('%s\n\t** The images do not all have same orientation and/or voxel sizes **',str);
    if verbose
        str = sprintf('%sThe function assumes that a voxel in one image  corresponds exactly\n',str);
        str = sprintf('%swith  the same voxel in another.   This is not a safe assumption if\n',str);
        str = sprintf('%sthe orientation information  in the headers or .mat files says that\n',str);
        str = sprintf('%sthe images are oriented differently. Please ensure that you process\n',str);
        str = sprintf('%sall data correctly. For example, you may have realigned the images,\n',str);
        str = sprintf('%sbut not actually resliced them to be in voxel-wise alignment.\n',str);
        str = sprintf('%sHere are the orientation matrices of the image volumes.   This list\n',str);
        str = sprintf('%scan be used to determine which file(s) are causing the problem.\n\n',str);
        for i=1:numel(V)
            str = sprintf('%s[%g %g %g %g; %g %g %g %g; %g %g %g %g]  %s\n',...
                str, V(i).mat(1:3,:)', V(i).fname);
        end
        str = sprintf('%s\n',str);
    end
    if ~nargout, error('The orientations etc must be identical for this procedure.'); end
end

if nargout<2
    disp(str);
end
