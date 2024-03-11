function MTsat = hmri_calc_MTsat(MTw, A, R1)
%hmri_calc_MTsat Calculate MTsat from MTw data.
% 
% Input:
%   MTw.data (array of MT-weighted signals)
%   MTw.fa   (nominal excitation flip angle, rad)
%   MTw.TR   (repetition time, s or ms)
%   MTw.B1   (array of B1 ratios: either 1 [classic MTsat] or actual fa / nominal fa [calibrated MTsat for high FA])
%   A        (unnormalised PD map, must not have been scaled relative to MTw.data)
%   R1       (R1 estimates, in reciprocal units of MTw.TR)
%
% Output:
%   MTsat (in percent units [p.u.])
%
% Examples:
%   Procedure from Helms (2008) correcting for B1 to second order within the small flip angle regime by setting B1 = 1:
%     A_for_MTsat = hmri_calc_A(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',1),...
%                               struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',1), true);
%     R1_for_MTsat = hmri_calc_R1(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',1),...
%                                 struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',1), true);
%     MTsat = hmri_calc_MTsat(struct('data',data_mtw,'fa',fa_mtw,'TR',tr_mtw,'B1',1), A_for_MTsat, R1_for_MTsat);
%
% References:
%   Helms, et al. (2008), "High-resolution maps of magnetization transfer with inherent correction for RF inhomogeneity 
%       and T1 relaxation obtained from 3D FLASH MRI". Magn. Reson. Med. https://doi.org/10.1002/mrm.21732
%   Tabelow, et al. (2019), "hMRI – A toolbox for quantitative MRI in neuroscience and clinical research". Neuroimage 
%       https://doi.org/10.1016/j.neuroimage.2019.01.029

if isempty(MTw.B1), MTw.B1=1; end

assert(all(size(MTw.data)==size(A)), 'hmri:inputArraySize','MTw.data and A must be the same size!')
assert(all(size(MTw.data)==size(R1)),'hmri:inputArraySize','MTw.data and R1 must be the same size!')
assert(MTw.TR>0,'hmri:TR','MTw.TR must be positive!')

MTsat = ( (A .* (MTw.B1*MTw.fa) - MTw.data) ./ (MTw.data+eps) .* R1*MTw.TR - (MTw.B1*MTw.fa).^2 / 2 ) * 100;

% Make data points with missing data NaN
nanmask=(MTw.data==0)|(MTw.data==0);
MTsat(nanmask)=NaN;

end