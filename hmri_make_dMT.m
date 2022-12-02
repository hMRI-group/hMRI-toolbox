function dMT = hmri_make_dMT(SPD,ST1,SMT,dSPD,dST1,dSMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT,small_angle_approximation)
% Calculate propagation of uncertainty for MT map
% (https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
% Taking the total differential from Eq. A9 in Tabelow et al., NI, 2019
% gives:
% S.Mohammadi 06.09.2019
%
% In:
% SPD           - PDw signal at TE=0
% ST1           - T1w signal at TE=0
% SMT           - MTw signal at TE=0
% dSPD          - residual of mono-exponential fit of PDw signal
% dST1          - residual of mono-exponential fit of T1w signal
% dSMT          - residual of mono-exponential fit of MTw signal
% alpha_PD      - flip angle of T1w signal
% alpha_MT      - flip angle of MTw signal
% TRPD          - repitition time of PDw signal
% TRT1          - repitition time of T1w signal
% TRMT          - repitition time of MTw signal
% small_angle_approximation - Switch to turn off small-angle approximation
% 
% Out:
% dMT           - error for MT in [a.u.]
%
% References:
%   Mohammadi et al. NeuroImage (2022), "Error quantification in 
%     multi-parameter mapping facilitates robust estimation and enhanced 
%     group level sensitivity." 
%     https://doi.org/10.1016/j.neuroimage.2022.119529

% We do not scale the flip angles by fT, because that would be inconsistent
% with what is used for MT calculation in the toolbox.

% If small angle approximation is not used for R1 and PD calculation, then
% we should not use it here either. Note that this does not affect
% alpha_MT, which is always assumed to be small enough to use the small
% angle approximation.
if ~small_angle_approximation
    alpha_PD=2*tan(alpha_PD/2);
    alpha_T1=2*tan(alpha_T1/2);
end

dMTdSPD = dMT_by_dS1(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT);
dMTdST1 = dMT_by_dS1(ST1,SPD,SMT,alpha_T1,alpha_PD,alpha_MT,TRT1,TRPD,TRMT);
dMTdSMT = dMT_by_dSMT(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT);

dMT = sqrt( dMTdSPD.^2 .* dSPD.^2 + dMTdST1.^2 .* dST1.^2 + dMTdSMT.^2 .* dSMT.^2);

end

function d = dMT_by_dS1(S1,S2,SMT,alpha1,alpha2,alphaMT,TR1,TR2,TRMT)
% Derivative of dual flip-angle A (PD) estimate with respect to first 
% weighted signal (S1). Because of symmetry in the R1 and A calculations, 
% the derivative with respect to the second weighted signal can be computed 
% by permuting labels.
%
% Can be derived using: 
%   syms S1 alpha1 S2 alpha2 SMT alphaMT
%   syms TR1 TR2 TRMT positive
%   A = hmri_calc_A(struct('data',S1,'fa',alpha1,'TR',TR1,'B1',1),struct('data',S2,'fa',alpha2,'TR',TR2,'B1',1),true);
%   R1 = hmri_calc_R1(struct('data',S1,'fa',alpha1,'TR',TR1,'B1',1),struct('data',S2,'fa',alpha2,'TR',TR2,'B1',1),true);
%   MT = (A*alphaMT/SMT-1)*R1*TRMT-alphaMT^2/2;
%   diff(MT,S1)

d = S2.*TRMT.*alpha1.*(TR2.*alpha1.^2 - TR1.*alpha2.^2).*(S2.*alphaMT - SMT.*alpha2)./(2.*SMT.*TR1.*TR2.*(S1.*alpha2 - S2.*alpha1).^2);

end

function d = dMT_by_dSMT(S1,S2,SMT,alpha1,alpha2,alphaMT,TR1,TR2,TRMT)
% Derivative of dual flip-angle A (PD) estimate with respect to MT 
% weighted signal (SMT).
%
% Can be derived using: 
%   syms S1 alpha1 S2 alpha2 SMT alphaMT
%   syms TR1 TR2 TRMT positive
%   A = hmri_calc_A(struct('data',S1,'fa',alpha1,'TR',TR1,'B1',1),struct('data',S2,'fa',alpha2,'TR',TR2,'B1',1),true);
%   R1 = hmri_calc_R1(struct('data',S1,'fa',alpha1,'TR',TR1,'B1',1),struct('data',S2,'fa',alpha2,'TR',TR2,'B1',1),true);
%   MT = (A*alphaMT/SMT-1)*R1*TRMT-alphaMT^2/2;
%   diff(MT,SMT)

d = S1.*S2.*TRMT.*alphaMT.*(TR2.*alpha1.^2 - TR1.*alpha2.^2)./(2*SMT.^2.*TR1.*TR2.*(S1.*alpha2 - S2.*alpha1));

end
