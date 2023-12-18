function [sts, str] = hmri_check_nifti_orientations(V, verbose)
% Check the dimensions and orientations of the images
% FORMAT [sts, str] = hmri_check_nifti_orientations(V [,verbose])
% V       - a struct array as returned by spm_vol
% verbose - [Default: true]
%
% sts     - status (true means OK)
% str     - string describing status, empty if OK
%
% When used without LHS, this function throws an error accordingly.
% If zero or one output arguments are requested, then str is output
% to the console.
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
    strdim = sprintf('\n\t** The images do not all have the same dimensions **\n');
    if verbose
        strdim = sprintf('%sThe function assumes that a voxel in one image corresponds with\n',strdim);
        strdim = sprintf('%sthe same  voxel in another.   This  is not a safe assumption if\n',strdim);
        strdim = sprintf('%sthe  image dimensions differ.   Please  ensure  that  you  have\n',strdim);
        strdim = sprintf('%sprocessed all the image data in the same way (eg. check spatial\n',strdim);
        strdim = sprintf('%snormalisation bounding-boxes, voxel-sizes etc).\n',strdim);
        strdim = sprintf('%sHere are the dimensions of the image volumes.  This list can be\n',strdim);
        strdim = sprintf('%sused to determine which file(s) are causing the problem.\n\n',strdim);
        for i=1:numel(V)
            strdim = sprintf('%s[%d %d %d]  %s\n', strdim, V(i).dim, V(i).fname);
        end
        strdim = sprintf('%s\n',strdim);
    end
    str = sprintf('%s%s',str,strdim);
    if nargout<2, disp(strdim); end
    if ~nargout, error('The dimensions must be identical for this procedure.'); end
end

matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4))
    sts = false;
    strori = sprintf('\n\t** The images do not all have same orientation and/or voxel sizes **\n');
    if verbose
        strori = sprintf('%sThe function assumes that a voxel in one image  corresponds exactly\n',strori);
        strori = sprintf('%swith  the same voxel in another.   This is not a safe assumption if\n',strori);
        strori = sprintf('%sthe orientation information  in the headers or .mat files says that\n',strori);
        strori = sprintf('%sthe images are oriented differently. Please ensure that you process\n',strori);
        strori = sprintf('%sall data correctly. For example, you may have realigned the images,\n',strori);
        strori = sprintf('%sbut not actually resliced them to be in voxel-wise alignment.\n',strori);
        strori = sprintf('%sHere are the orientation matrices of the image volumes.   This list\n',strori);
        strori = sprintf('%scan be used to determine which file(s) are causing the problem.\n\n',strori);
        for i=1:numel(V)
            strori = sprintf('%s[%g %g %g %g; %g %g %g %g; %g %g %g %g]  %s\n',...
                strori, V(i).mat(1:3,:)', V(i).fname);
        end
        strori = sprintf('%s\n',strori);
    end
    str = sprintf('%s%s',str,strori);
    if nargout<2, disp(strori); end
    if ~nargout, error('The orientations etc must be identical for this procedure.'); end
end
