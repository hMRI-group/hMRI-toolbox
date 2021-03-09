% ========================================================================
% function hmri_read_vols(V,VG,p,interp)
%
% read image volume
%
% Input:
%   V      - structure containing image volume information of ith image
%   VG     - structure containing image volume information of target image
%   p      - z position
%   interp - interpolation
%   x      - output of spm_coreg for when an additional coregistration step
%            is required between V and VG, e.g. taking echoes from native
%            space but putting the output in an echo-averaged space that
%            has already been coregisterd to V_pdw(1) giving x.
%
% ========================================================================
% S.Mohammadi 18/10/2019

function dataOut = hmri_read_vols(V,VG,p,interp, x)

if ~exist('x', 'var')
    x = zeros(1,6);
end
dm = VG.dim;
% M = spm_matrix([0 0 p 0 0 0 1 1 1]);
% M1 = V.mat\VG.mat*M;
% dataOut = spm_slice_vol(V,M1,dm(1:2),interp);

M = inv(V.mat)*spm_matrix(x)*VG.mat*spm_matrix([0 0 p 0 0 0]);
dataOut = spm_slice_vol(V,M,dm(1:2),interp);

end