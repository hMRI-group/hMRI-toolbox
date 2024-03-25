function fn_out = hmri_proc_zero2nan(fn_in)
% This function turns 0's into NaN's in the input image(s) when their data 
% format allows it, i.e. for images in float (32 or 64) format. The 0's 
% naturally appear when masking out images with a binary mask, i.e. when 
% mulitplying voxels in an image with 1 or 0, e.g. when (tissue-weighted) 
% smoothing spatially normalized quantitative maps. 
% 
% The issue is that in SPM, when analyzing float images, their 0's are NOT 
% implictly masked out during a GLM analysis and voxels with a 0 are thus 
% considered a genuine signal (if not masked out explicitly).
% This could bias the stats at those voxels and lead to spurious 
% (possibly significant) effects; e.g take a voxel being 0 moreoften in one
% group than the other, then group mean values will likely be different.
% 
% This routine will thus simply rewrite the input image(s) with NaN instead 
% of 0's. The routine will check the data format and, if it's in float,
% turn all the 0's into NaN's. Other values are left untouched.
% NOTE: the image values are OVERWRITTEN and output list of filenames
% will be the same as the input.
% 
% 
% FORMAT
%   fn_out = hmri_proc_zero2nan(fn_in)
% 
% INPUT 
%   fn_in   : list (char or cell array) of image filenames
% 
% OUTPUT 
%   fn_out  : list (char or cell array) of image filenames, same as fn_in
% 
% LIMITATION
% This only works on 3D images.
% A list of frames, i.e. 3D volumes, from a 4D image can be entered though.
% 
% NOTE
% The routine can also be used to fix a bunch of previously smoothed
% quantitative maps, aka. wap1/2* images, to ensure that 0's are turned
% into NaN's and those voxels implicitly masked out when entering them into
% a GLM in SPM
% Example:
% 
%   fn_waps = spm_select('FPList','path_to_my_data_folder','^wap_.*\.nii$')
%   hmri_proc_zero2nan(fn_waps)
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Prepare input/output
% Images are overwritten so same output file names
fn_out = fn_in;

% Turn a cell array into char array
if iscell(fn_in)
    fn_in = char(fn_in);
end

% Number of images to deal with
Nimg = size(fn_in,1);

%% Deal with all the images.
% Since current computer have large memories, we can proceed directly by
% using the file_array object to replace 0's with NaN's, when needed.
for ii=1:Nimg
    V_ii = spm_vol(deblank(fn_in(ii,:))); % memory map each image
    if spm_type(V_ii.dt(1),'nanrep')
        % Do the conversion if it has NaN representation
        dd = V_ii.private.dat(:,:,:); % read in the whole 3D image
        dd(dd(:)==0) = NaN; % Turn 0's into NaN in data array
        V_ii.private.dat(:,:,:) = dd; % write back data into file
    end
end

end
