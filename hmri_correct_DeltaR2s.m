function DeltaR2s = hmri_correct_DeltaR2s(DeltaR2s, f_T)
% Correct DeltaR2s for local flip angle inhomogeneity. This is not done in
% hmri_calc_R2s, as often the f_T map is not available when R2* is being
% fitted.
%
% As for now, only a linear dependence on flip angle is fitted by
% hmri_calc_R2s, we just divide out the inhomogeneity. If other methods
% are implemented then different corrections will also need to be added
% here.
%
% Reference:
%   Milotta et al. Magn. Reson. Med. (2023), "Mitigating the impact of
%     flip angle and orientation dependence in single compartment R2*
%     estimates via 2-pool modeling."
%     https://doi.org/10.1002/mrm.29428

DeltaR2s = DeltaR2s./f_T;

end