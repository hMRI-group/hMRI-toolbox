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
%
% ========================================================================
% S.Mohammadi 18/10/2019

function dataOut = hmri_read_vols(V,VG,p,interp)

dm = VG.dim;
M = spm_matrix([0 0 p 0 0 0 1 1 1]);
M1 = V.mat\VG.mat*M;
dataOut = spm_slice_vol(V,M1,dm(1:2),interp);

end