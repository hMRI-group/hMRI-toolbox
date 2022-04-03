function [dPD,AdPD] = hmri_make_dPD(SPD,ST1,dSPD,dST1,alpha_PD,alpha_T1,TRPD,TRT1,A,f_T,threshall,small_angle_approximation)
% Calculate propagation of uncertainty for PD map
% (https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
% Here the code for generating the derivatives of Eq. A7 in Tabelow et al., NI, 2019:
% syms PD SPD ST1 alpha_PD alpha_T1 TRPD TRT1
% PD = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) ST1.*SPD.*(TRT1*alpha_PD/alpha_T1-TRPD*alpha_T1/alpha_PD)./(TRPD*SPD*alpha_PD-TRT1*ST1*alpha_T1);
% diff(PD,SPD)
%                       / TRPD alpha_T1   TRT1 alpha_PD \
% SPD ST1 TRPD alpha_PD | ------------- - ------------- |
%                       \    alpha_PD        alpha_T1   /
% -------------------------------------------------------
%                                                2
%         (SPD TRPD alpha_PD - ST1 TRT1 alpha_T1)
%
%          / TRPD alpha_T1   TRT1 alpha_PD \
%      ST1 | ------------- - ------------- |
%          \    alpha_PD        alpha_T1   /
%    - -------------------------------------
%      SPD TRPD alpha_PD - ST1 TRT1 alpha_T1
% diff(PD,ST1)
%       / TRPD alpha_T1   TRT1 alpha_PD \
%   SPD | ------------- - ------------- |
%       \    alpha_PD        alpha_T1   /
% - -------------------------------------
%   SPD TRPD alpha_PD - ST1 TRT1 alpha_T1
%
%                            / TRPD alpha_T1   TRT1 alpha_PD \
%      SPD ST1 TRT1 alpha_T1 | ------------- - ------------- |
%                            \    alpha_PD        alpha_T1   /
%    - -------------------------------------------------------
%                                                     2
%              (SPD TRPD alpha_PD - ST1 TRT1 alpha_T1)
%
% % S.Mohammadi 06.09.2019
%
% In:
% SPD           - PDw signal at TE=0
% ST1           - T1w signal at TE=0
% dSPD          - residual of mono-exponential fit of PDw signal
% dST1          - residual of mono-exponential fit of T1w signal
% alpha_PD      - flip angle of PDw signal
% alpha_T1w     - flip angle of T1w signal
% TRPD          - repitition time of PDw signal
% TRT1          - repitition time of T1w signal
% A_forMT       - apparent proton density in arbritrary units
% f_T           - map of transmit field
% VG            - target structure of
%
% Out:
% dPD           - error for A in [a.u.]
% AdPD          - error map for A in [a.u.]

if(~isempty(f_T))
    alpha_PD = alpha_PD.*f_T;
    alpha_T1 = alpha_T1.*f_T;
    if ~small_angle_approximation
        alpha_PD=2*tau(alpha_PD/2);
        alpha_T1=2*tau(alpha_T1/2);
    end
end

% dPD calculation is symmetric with respect to the two weighted contrasts
dPD = sqrt( dPD_by_dS1(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dSPD.^2 ...
    +dPD_by_dS1(ST1,SPD,alpha_T1,alpha_PD,TRT1,TRPD).^2.*dST1.^2 );

AdPD     = zeros(size(SPD));
tmp1    = dPD;

% TODO: input argument "A" is unused
% TODO: handle thresholding outside of this function
A = max(min(tmp1,threshall.A),-threshall.A);
AdPD(A>threshall.dPD) = tmp1(A>threshall.dPD);

end

function d = dPD_by_dS1(S1,S2,alpha1,alpha2,TR1,TR2)
d = S1.*S2.*TR1.*alpha1.*(TR1.*alpha2./alpha1 - TR2.*alpha1./alpha2) ./ (S1.*TR1.*alpha1 - S2.*TR2.*alpha2).^2 ...
    - S2.*(TR1.*alpha2./alpha1 - TR2.*alpha1./alpha2)./(S1.*TR1.*alpha1 - S2.*TR2.*alpha2);
end