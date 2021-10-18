function [R2s,extrapolated]=hmri_calc_R2s(weighted_data,method)
% R2* estimation using the ESTATICS (ESTimating the Apparent Transverse 
% relaxation time from Images with different ContrastS) model
%
% array of structures (one per contrast) input in the form:
%   weighted(contrast).data (NvoxelsX x NvoxelsY x ... x Nechoes)
%   weighted(contrast).TE   (1 x Nechoes)
% Voxels must correspond between the weightings (i.e. the images should 
% have been resliced to the same space), but the sampled TEs may be 
% different.
% Nechoes must be at least 2 for each weighting.
% 
% outputs:
%   R2s (NvoxelsX x NvoxelsY x ...) contains the voxelwise-estimated 
%       common R2* of the weightings.
%   extrapolated contains size(weighted_data) extrapolated to TE=0 in 
%       the same order as the input (e.g. matching contrast order).

assert(isstruct(weighted_data(1)),'hmri:structError',['inputs must be structs; see help ' mfilename])

switch lower(method)
    case {'ols','robust','wls1','wls','wls3'}
        [R2s,extrapolated]=hmri_calc_R2sLL(weighted_data,method);        
    case {'arlo','darlo','nlls_ols','nlls_wls1'}
        [R2s,extrapolated]=hmri_calc_R2sNL(weighted_data,method);
    otherwise
        error('fitting method not available!')
end
