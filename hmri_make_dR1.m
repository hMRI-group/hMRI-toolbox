function [dR1,AdR1] = hmri_make_dR1(SPD,ST1,dSPD,dST1,alpha_PD,alpha_T1,TRPD,TRT1,f_T,R1,VG,threshall)
% Calculate propagation of uncertainty for R1 map
% (https://en.wikipedia.org/wiki/Propagation_of_uncertainty). 
% Here the code for generating the derivatives of Eq. A6 in Tabelow et al., NI, 2019:
% syms R1 SPD ST1 alpha_PD alpha_T1 TRPD TRT1
% R1 = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) 0.5*(SPD.* alpha_PD/TRPD - ST1.* alpha_T1/TRT1)./(ST1./alpha_T1 - SPD./alpha_PD);
% diff(R1,SPD)=
%    SPD alpha_PD   ST1 alpha_T1
%    ------------ - ------------
%       2 TRPD         2 TRT1                    alpha_PD
% --------------------------------- - ------------------------------
%          /    SPD        ST1   \2        /    SPD        ST1   \
% alpha_PD | -------- - -------- |    TRPD | -------- - -------- | 2
%          \ alpha_PD   alpha_T1 /         \ alpha_PD   alpha_T1 /
% diff(R1,ST1)=
%                                     SPD alpha_PD   ST1 alpha_T1
%                                     ------------ - ------------
%            alpha_T1                    2 TRPD         2 TRT1
% ------------------------------ - ---------------------------------
%      /    SPD        ST1   \              /    SPD        ST1   \2
% TRT1 | -------- - -------- | 2   alpha_T1 | -------- - -------- |
%      \ alpha_PD   alpha_T1 /              \ alpha_PD   alpha_T1 /
% S.Mohammadi 06.09.2019
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
% f_T           - map of transmit field 
% R1            - longitudinal relaxation rate parameter map
% VG            - target structure of 
% 
% Out:
% dR1           - error for R1 in [1/ms] 
% AdR1          - error map for R1 in [1/s]

dm = VG.dim;
if(~isempty(f_T))
    alpha_PD = alpha_PD.*f_T;
    alpha_T1 = alpha_T1.*f_T;
end  

dR1dSPD = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) ((SPD.*alpha_PD)./(2*TRPD) - (ST1.*alpha_T1)./(2*TRT1))./(alpha_PD.*(SPD./alpha_PD - ST1./alpha_T1).^2) - alpha_PD./(2*TRPD*(SPD./alpha_PD - ST1./alpha_T1));
dR1dST1 = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) alpha_T1./(2*TRT1.*(SPD./alpha_PD - ST1./alpha_T1)) - ((SPD.*alpha_PD)./(2*TRPD) - (ST1.*alpha_T1)./(2*TRT1))./(alpha_T1.*(SPD./alpha_PD - ST1./alpha_T1).^2);

dR1 = sqrt(dR1dSPD(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dSPD.^2+dR1dST1(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dST1.^2);

tmp     = zeros(dm(1:2));
tmp1    = dR1*1e3;
tmp1 = max(min(tmp1,threshall.R1),-threshall.R1);
tmp(R1>threshall.dR1)     = tmp1(R1>threshall.dR1);
AdR1 = tmp;