% ========================================================================
% function x = hmri_coreg(P_ref, P_src, coreg_flags)
%
% coregistration of volumes
%
% Input:
%   P_ref  - input to spm_vol that wil be reference for co-registration
%   P_src  - input to spm_vol that wil be co-registered (source)
%   coreg_flags - flags for co-registration, e.g. interpolation
%
% Output:
%   x      - output of spm_coreg.
%
% ========================================================================

function x = hmri_coreg(P_ref, P_src, coreg_flags)

assert(size(P_src, 1) <= 2, ...
    'Expecting a maximum of two volumes in P_src to be co-registered to P_ref');

VG = spm_vol(P_ref);
VF = spm_vol(P_src(1,:));
x = spm_coreg(VG,VF,coreg_flags);
M  = inv(spm_matrix(x));
MM = spm_get_space(deblank(VF.fname));
spm_get_space(deblank(deblank(VF.fname)), M*MM);

if size(P_src, 1) == 2
    % This should be a B1 map, apply transform to it too
    VF2 = spm_vol(P_src(2,:));
    M  = inv(spm_matrix(x));
    MM = spm_get_space(deblank(VF2.fname));
    spm_get_space(deblank(deblank(VF2.fname)), M*MM);
end

    
end