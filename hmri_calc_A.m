function A=hmri_calc_A(PDw,T1w,small_angle_approx)
%hmri_calc_A Calculate A (an unnormalised PD map) from PDw and T1w data.
% 
% Input:
%   PDw and T1w are structs with the fields:
%     {PDw,T1w}.data (array of signals)
%     {PDw,T1w}.fa   (nominal flip angle, rad)
%     {PDw,T1w}.TR   (repetition time, s or ms)
%     {PDw,T1w}.B1   (array of B1 ratios: actual fa / nominal fa)
%
%   small_angle_approx (bool which determines whether to use the small
%     angle approximation [true] or to use a more exact formula [false])
%
% Output:
%   A (in arbitrary units)
%
% Examples:
%   Estimate A using the small angle approximation:
%       A = hmri_calc_A(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',b1map),...
%            struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',b1map), true);
%
%   Estimate A without the small angle approximation:
%       A = hmri_calc_A(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',b1map),...
%            struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',b1map), false);
%
% References:
%   Helms et al. Magn. Reson. Med. (2008), "Quantitative FLASH MRI at 3T 
%       using a rational approximation of the Ernst equation",
%       https://doi.org/10.1002/mrm.21542
%
%   If you use small_angle_approx=false:
%       Edwards et al.  Magn. Reson. Mater. Phy. (2021), "Rational 
%           approximation of the Ernst equation for dual angle R1 mapping 
%           revisited: beyond the small flip-angle assumption" in Book of 
%           Abstracts ESMRMB 2021, 
%           https://doi.org/10.1007/s10334-021-00947-8

if isempty(PDw.B1), PDw.B1=1; end
if isempty(T1w.B1), T1w.B1=1; end

if ~small_angle_approx
    PDw.t=2*tan(PDw.B1.*PDw.fa/2);
    T1w.t=2*tan(T1w.B1.*T1w.fa/2);
else
    PDw.t=PDw.B1.*PDw.fa;
    T1w.t=T1w.B1.*T1w.fa;
end

assert(all(size(PDw.data)==size(T1w.data)),'hmri:inputArraySize','PDw.data and T1w.data must be the same size!')

% Sanity check TR
% Use "isAlways" so that hmri_calc_A.m can be called using symbolic variables
assert(isAlways(PDw.TR>0),'hmri:TR','PDw.TR must be positive!')
assert(isAlways(T1w.TR>0),'hmri:TR','T1w.TR must be positive!')

A=T1w.data.*PDw.data.*( T1w.TR.*PDw.t./T1w.t - PDw.TR.*T1w.t./PDw.t )...
    ./ ( PDw.data.*T1w.TR.*PDw.t - T1w.data.*PDw.TR.*T1w.t );

% Make data points with missing data NaN
nanmask=(T1w.data==0)|(PDw.data==0)|(T1w.t==0)|(PDw.t==0);
A(nanmask)=NaN;

end
