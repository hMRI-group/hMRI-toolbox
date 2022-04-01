function [dPD,Atmp] = hmri_make_dPD(SPD,ST1,dSPD,dST1,alpha_PD,alpha_T1,TRPD,TRT1,A,VG,f_T,threshall)
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

dm = VG.dim;
% TODO: include version without small angle approximation for PD
% calculation
if(~isempty(f_T))
    alpha_PD = alpha_PD.*f_T;
    alpha_T1 = alpha_T1.*f_T;
end  

dPDSPD = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) (SPD.*ST1.*TRPD.*alpha_PD.*((TRPD.*alpha_T1)./alpha_PD - (TRT1.*alpha_PD)./alpha_T1))./(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1).^2 - (ST1.*((TRPD.*alpha_T1)./alpha_PD - (TRT1.*alpha_PD)./alpha_T1))./(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1);

dPDSR1 = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) - (SPD.*((TRPD.*alpha_T1)./alpha_PD - (TRT1*alpha_PD)./alpha_T1))./(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1) - (SPD.*ST1.*TRT1.*alpha_T1.*((TRPD*alpha_T1)./alpha_PD - (TRT1.*alpha_PD)./alpha_T1))./(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1).^2;
dPD = sqrt(dPDSPD(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dSPD.^2+dPDSR1(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dST1.^2);
Atmp     = zeros(dm(1:2));
tmp1    = dPD;

% TODO: input argument "A" is unused
A = max(min(tmp1,threshall.A),-threshall.A);
Atmp(A>threshall.dPD)     = tmp1(A>threshall.dPD);
