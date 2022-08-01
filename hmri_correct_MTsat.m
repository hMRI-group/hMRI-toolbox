function MTsat=hmri_correct_MTsat(MTsat,B1)
%hmri_correct_MTsat Correct MTsat using B1 and a heuristic model.
% 
% Input:
%   PDw and T1w are structs with the fields:
%     MTsat (array of MTsat estimates)
%     B1    (array of B1 ratios: actual fa / nominal fa)
%
% Output:
%   MTsat (in percent units [p.u.])
%
% Examples:
% References:

MTsat = MTsat .* (1 - 0.4) ./ (1 - 0.4 * B1);

% Make data points with missing data NaN
nanmask=(MTsat.data==0)|(MTsat.data==0);
MTsat(nanmask)=NaN;

end
