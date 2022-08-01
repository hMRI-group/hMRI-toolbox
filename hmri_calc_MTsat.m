function MTsat=hmri_calc_MTsat(MTw,A,R1)
%hmri_calc_MTsat Calculate MTsat from MTw data.
% 
% Input:
%     MTw.data (array of MT-weighted signals)
%     MTw.fa   (nominal flip angle, rad)
%     MTw.TR   (repetition time, s or ms)
%     MTw.B1   (array of B1 ratios: either 1 [classic MTsat] or actual fa / nominal fa [calibrated MTsat for high FA])
%     A        (unnormalised PD map, must not have been scaled relative to MTw.data)
%     R1       (R1 estimates, in reciprocal units of MTw.TR)
%
% Output:
%   MTsat (in percent units [p.u.])
%
% Examples:
% References:

if isempty(MTw.B1), MTw.B1=1; end

assert(all(size(MTw.data)==size(A.data)),'hmri:inputArraySize','MTw.data and A must be the same size!')
assert(all(size(MTw.data)==size(R1.data)),'hmri:inputArraySize','MTw.data and R1 must be the same size!')
assert(MTw.TR>0,'hmri:TR','MTw.TR must be positive!')

MTsat = ( (A .* (MTw.B1*MTw.fa) - MTw.data) ./ (MTw.data+eps) .* R1*MTw.TR - (MTw.B1*MTw.fa).^2 / 2 ) * 100;

% Make data points with missing data NaN
nanmask=(MTw.data==0)|(MTw.data==0);
MTsat(nanmask)=NaN;

end
