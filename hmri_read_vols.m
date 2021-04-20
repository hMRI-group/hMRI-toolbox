% ========================================================================
% function hmri_read_vols(V,VG,p,interp,x)
%
% read image volume
%
% Input:
%   V      - structure containing image volume information of ith image
%   VG     - structure containing image volume information of target image
%   p      - z position
%   interp - interpolation for spm_vol. Values between -127 and 127
%   x      - optional argument that is the output of spm_coreg.  For use 
%            when an additional coregistration step
%            is required between V and VG, e.g. taking echoes from native
%            space but putting the output in an echo-averaged space that
%            has already been coregisterd to V_pdw(1) giving x.
%
% ========================================================================
% S.Mohammadi 18/10/2019

function dataOut = hmri_read_vols(V,VG,p,interp,x)

assert(isstruct(V),'hmri:structError',['Input V must be struct; see help ' mfilename])
assert(isstruct(VG),'hmri:structError',['Input VG must be struct; see help ' mfilename])
assert(isscalar(interp) && abs(interp) < 128,'hmri:inputError','Invalid interpolation setting; see help spm_slice_vol')
assert(isnumeric(p),'hmri:typeError',['z position input must be numerical; see help ' mfilename])

if ~exist('x', 'var') || isempty(x)
    x = zeros(1,6);
else
    assert(isvector(x),'hmri:typeError','x must be vector output from spm_coreg; see help spm_coreg')
    x = x(:)';
end
    
dm = VG.dim;
M = inv(V.mat)*spm_matrix(x)*VG.mat*spm_matrix([0 0 p 0 0 0]);
dataOut = spm_slice_vol(V,M,dm(1:2),interp);

end