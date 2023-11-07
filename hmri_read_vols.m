% ========================================================================
% function hmri_read_vols(V,VG,p,interp,x)
%
% read image volume
%
% Input:
%   V      - structure from spm_vol containing image volume information 
%            of image to be read
%   VG     - structure from spm_vol or SPM nifti object containing image 
%            volume information of image defining the target space
%   p      - z position
%   interp - interpolation value for spm_vol. Values between -127 and 127
%   x      - optional argument that is the output of spm_coreg. For use 
%            when an additional transformation is required between V and 
%            VG, e.g. when V is defined in native space but the transform
%            from V space to VG space has been previously estimated using 
%            another image V' (e.g. the average over echoes) in the same 
%            space as V.
%
% ========================================================================
% S.Mohammadi 18/10/2019

function dataOut = hmri_read_vols(V,VG,p,interp,x)

assert(isstruct(V),'hmri:structError',['Input V must be struct from spm_vol; see help ' mfilename])
assert(isstruct(VG)||isa(VG,'nifti'),'hmri:structError',['Input VG must be struct from spm_vol or SPM nifti object; see help ' mfilename])
assert(isscalar(interp) && abs(interp) < 128,'hmri:inputError','Invalid interpolation setting; see help spm_slice_vol')
assert(isnumeric(p),'hmri:typeError',['z position input must be numerical; see help ' mfilename])

if ~exist('x', 'var') || isempty(x)
    x = zeros(1,6);
else
    assert(isvector(x),'hmri:typeError','x must be vector output from spm_coreg; see help spm_coreg')
    x = x(:)';
end
    
if isa(VG,'nifti')
    dm = VG.dat.dim;
elseif isstruct(VG)
    dm = VG.dim;
end

M = inv(V.mat)*spm_matrix(x)*VG.mat*spm_matrix([0 0 p 0 0 0]);
dataOut = spm_slice_vol(V,M,dm(1:2),interp);

end
