% ========================================================================
% function x = hmri_coreg(P_ref, P_src, flags)
%
% coregistration of volumes
%
% Input:
%   P_ref - input to spm_vol that will be reference for co-registration
%   P_src - input to spm_vol that will be co-registered (source); The 
%           first entry will always be used to determine the
%           transformation. If P_src contains multiple entries (e.g. 
%           anatomical (first) + B1 map) then the header of each will be 
%           updated.
%   flags - co-registration flags, e.g. interpolation, and brain masking
%           options
%
% Output:
%   x     - output of spm_coreg.
%
% ========================================================================

function x = hmri_coreg(P_ref, P_src, flags)

assert(~isempty(P_ref), 'P_ref must not be empty');
assert(~isempty(P_src), 'P_src must not be empty');

% Transformation estimation:
VG = spm_vol(P_ref);
VF = spm_vol(P_src(1,:));
if isfield(flags, 'mask_options') && flags.mask_options.domask
    [~,Vmask] = hmri_create_pm_brain_mask(VF, flags.mask_options.flags);
    if spm_type(VF.dt(1), 'nanrep')
        % hack for replacing 0 in mask with nan
        f = 'i1.*(i2./i2)';
    else
        % zero-out voxels outside the mask
        f = 'i1.*i2';
    end
    VF = spm_imcalc([VF, Vmask], spm_file(VF.fname, 'suffix', '_masked'), f, struct('dtype', VF.dt(1)));
else % use unmasked source image for registration
    VF = spm_vol(P_src(1,:));
end
x = spm_coreg(VG,VF,flags);

% Application of transformation to header (no reslicing):
M  = pinv(spm_matrix(x));
for ind = 1 : size(P_src, 1)
    VF = spm_vol(P_src(ind,:));
    MM = spm_get_space(deblank(VF.fname));
    spm_get_space(deblank(deblank(VF.fname)), M*MM);
end

end
