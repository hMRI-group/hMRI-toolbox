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
%   flags - co-registration flags, e.g. interpolation; can be empty.
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
x = spm_coreg(VG,VF,flags);

% Application of transformation to header (no reslicing):
M  = pinv(spm_matrix(x));
for ind = 1 : size(P_src, 1)
    VF = spm_vol(P_src(ind,:));
    MM = spm_get_space(deblank(VF.fname));
    spm_get_space(deblank(deblank(VF.fname)), M*MM);
end

end