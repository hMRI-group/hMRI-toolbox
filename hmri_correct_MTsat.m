function MTsat = hmri_correct_MTsat(MTsat, B1, model, C)
%hmri_correct_MTsat Correct MTsat using B1 and a heuristic model.
% 
% Input:
%   MTsat (array of MTsat estimates in percent units [p.u.])
%   B1    (array of B1 ratios: actual fa / nominal fa)
%   model (either 'lipp' or 'helms' model)
%   C     (parameter of heuristic model)
%
% Output:
%   MTsat (in p.u.)
%
% Examples:
%   If using the Helms, et al. (2021) model, then the input MTsat needs to have been computed with B1=1:
%     A_for_MTsat  = hmri_calc_A( struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',1), ...
%                                 struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',1), true);
%     R1_for_MTsat = hmri_calc_R1(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',1), ...
%                                 struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',1), true);
%     MTsat = hmri_calc_MTsat(struct('data',data_mtw,'fa',fa_mtw,'TR',tr_mtw,'B1',1), A_for_MTsat, R1_for_MTsat);
%   as the B1 correction then happens here:
%     MTsat_corrected = hmri_correct_MTsat(MTsat, B1, 'helms', 0.4);
%
%   If using the Lipp, et al. (2023) model, then the input MTsat needs to have been computed with the real B1. Note
%   that the T1w excitation flip angle was large in Lipp, et al. (2023), so the small angle approximation was not
%   used when computing A and R1.
%     A_for_MTsat  = hmri_calc_A( struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',B1), ...
%                                 struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',B1), false);
%     R1_for_MTsat = hmri_calc_R1(struct('data',data_pdw,'fa',fa_pdw,'TR',tr_pdw,'B1',B1), ...
%                                 struct('data',data_t1w,'fa',fa_t1w,'TR',tr_t1w,'B1',B1), false);
%     MTsat = hmri_calc_MTsat(struct('data',data_mtw,'fa',fa_mtw,'TR',tr_mtw,'B1',B1), A_for_MTsat, R1_for_MTsat);
%     MTsat_corrected = hmri_correct_MTsat(MTsat, B1, 'lipp', 1.2);
%
% References:
%   Helms, et al. (2021) "Correction of FLASH-based MT saturation in human brain for residual bias of B1-inhomogeneity at 3T"
%       arXiv preprint. https://doi.org/10.48550/arXiv.2104.14878
%   Lipp, et al. (2023). "B1+-correction of magnetization transfer saturation maps optimized for 7T postmortem MRI of the brain".
%       Magn. Reson. Med. https://doi.org/10.1002/mrm.29524

switch model
    case 'helms'
        MTsat = MTsat .* (1 - C) ./ (1 - C * B1);
    case 'lipp'
        MTsat = MTsat ./ (1 + C * (B1 - 1));
    otherwise
        error('unknown MTsat correction model ''%s''. Allowed models are ''helms'' and ''lipp''',model)
end

end