function fn_out = hmri_proc_zero2nan(fn_in)
% This is a "fix function" to deal with previously generate wap* images,
% i.e. tissue-weighted smoothed qMR images.
% Because the wap* images are in float, their 0's are NOT masked out during
% an SPM analysis and voxels with a 0 are thus considered a genuine signal.
% This could bias the stats at those voxels, if not masked explicitly, and
% lead to spurious (significant) results (think of a voxel being 0 more
% often in a group than another, the group means will likely be different).
% 
% The aim of this routine is simply to rewrite the wap* images with NaN
% instead of 0's. Note that this is already done explicitly with the 
% "latest" version of the smoothing function.
% 
% The routine will check the data format and if it's in float make sure
% that 0's are turned into NaN's. Other values are left untouched,
% therefore the image values are ovewritten and output list of filenames
% will be the same as the input.
% 
% FORMAT
%   fn_out = hmri_proc_zero2nan(fn_in)
% 
% INPUT 
%   fn_in   : list (char or cell array) of image filenames
% 
% OUTPUT 
%   fn_out  : list (char or cell array) of image filenames
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Prepare input/output
fn_out = fn_in;

if iscell(fn_in)
    fn_in = char(fn_in);
end
Nimg = size(fn_in,1);

%% Deal with all the images.
% Since current computer have large memories, we can proceed directly by
% using the file_array object to replace 0's with NaN's, when needed.
for ii=1:Nimg
    V_ii = spm_vol(deblank(fn_in(ii,:)));
    if spm_type(V_ii.dt(1),'nanrep')
        % Do the conversion if it has NaN representation
        dd = V_ii.private.dat(:,:,:);
        dd(dd(:)==0) = NaN;
        V_ii.private.dat(:,:,:) = dd;
    end
end

end
