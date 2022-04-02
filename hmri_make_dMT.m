function [dMT,AdMT] = hmri_make_dMT(SPD,ST1,SMT,dSPD,dST1,dSMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT,threshall,small_angle_approximation)
% old version: hmri_make_dMT(SMT,PD,R1,dSMT,dPD,dR1,alpha_MT,TRMT)
% old version: function [dMT,Atmp] = hmri_make_dMT(SMT,A_forMT,R1,MT,dSMT,dPD,dR1,alpha_MT,TRMT,VG)

% Calculate propagation of uncertainty for MT map
% (https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
% Taking the total differential from Eq. A9 in Tabelow et al., NI, 2019
% gives:
% S.Mohammadi 06.09.2019
%
% In:
% SMT           - MTw signal at TE=0
% A_forMT       - apparent proton density in arbritrary units
% R1            - relaxation rate in 1/ms(?)
% MT            - magnetization transfer saturation rate
% dSMT          - residual of mono-exponential fit of MTw signal
% dSPD          - residual of mono-exponential fit of PDw signal
% dST1          - residual of mono-exponential fit of T1w signal
% alpha_MT      - flip angle of MTw signal
% TRMT          - repitition time of MTDw signal
% VG            - target structure of
%
% Out:
% dMT           - error for MT in [a.u.]
% Atmp          - error map for MT in [a.u.]
% syms SPD ST1 alpha_PD alpha_T1 alpha_MT TRPD TRT1 TRMT SMT
% R1 = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) 0.5*(SPD.* alpha_PD/TRPD - ST1.* alpha_T1/TRT1)./(ST1./alpha_T1 - SPD./alpha_PD);
% Astar = @(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1) ST1.*SPD.*(TRT1*alpha_PD/alpha_T1-TRPD*alpha_T1/alpha_PD)./(TRPD*SPD*alpha_PD-TRT1*ST1*alpha_T1);
% MT = @(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT) (Astar(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1)*alpha_MT/SMT - 1)*R1(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1)*TRMT - alpha_MT^2/2;
%
% pretty(simplify(diff(MT,SPD)))
%
% (ST1 TRMT alpha_PD (TRPD alpha_T1  - TRT1 alpha_PD ) (alpha_MT SPD  ST1 TRPD  alpha_T1  - alpha_MT SPD  ST1 TRPD TRT1
%
%            2               2         2         2          2     2                  2                   2     2
%    alpha_PD  - alpha_MT SPD  ST1 TRT1  alpha_T1  + SMT SPD  TRPD  alpha_T1 alpha_PD  + alpha_MT SPD ST1  TRT1  alpha_T1
%
%                                                 2                        3                   2
%    alpha_PD 2 - 2 SMT SPD ST1 TRPD TRT1 alpha_T1  alpha_PD - alpha_MT ST1  TRPD TRT1 alpha_T1
%
%             2     2         3                                                           2                              2
%    + SMT ST1  TRT1  alpha_T1 ))/(2 SMT TRPD TRT1 (SPD TRPD alpha_PD - ST1 TRT1 alpha_T1)  (SPD alpha_T1 - ST1 alpha_PD) )
%
%
% pretty(simplify(diff(MT,ST1)))
%
% -(SPD TRMT alpha_T1 (TRPD alpha_T1  - TRT1 alpha_PD ) (- alpha_MT SPD  TRPD TRT1 alpha_PD
%
%                  2         2                              2     2         3                   2     2         2
%    + alpha_MT SPD  ST1 TRPD  alpha_T1 alpha_PD 2 + SMT SPD  TRPD  alpha_PD  - alpha_MT SPD ST1  TRPD  alpha_PD
%
%                      2                   2                   2     2         2
%    - alpha_MT SPD ST1  TRPD TRT1 alpha_T1  + alpha_MT SPD ST1  TRT1  alpha_PD  - 2 SMT SPD ST1 TRPD TRT1 alpha_T1
%
%            2          2     2         2                                                                    2
%    alpha_PD  + SMT ST1  TRT1  alpha_T1  alpha_PD))/(2 SMT TRPD TRT1 (SPD TRPD alpha_PD - ST1 TRT1 alpha_T1)
%
%                                 2
%    (SPD alpha_T1 - ST1 alpha_PD) )
%
%
%  pretty(simplify(diff(MT,SMT)))
%
%   SPD ST1 TRMT alpha_MT (SPD TRT1 alpha_PD - ST1 TRPD alpha_T1) (TRPD alpha_T1  - TRT1 alpha_PD )
% - -----------------------------------------------------------------------------------------------
%             2
%        2 SMT  TRPD TRT1 (SPD TRPD alpha_PD - ST1 TRT1 alpha_T1) (SPD alpha_T1 - ST1 alpha_PD)

% We do not scale the flip angles by fT, because that would be inconsistent
% with what is used for MT calculation in the toolbox.

% If small angle approximation is not used for R1 and PD calculation, then
% we should not use it here either. Note that this does not affect
% alpha_MT, which is always assumed to be small enough to use the small
% angle approximation.
if ~small_angle_approximation
    alpha_PD=2*tan(alpha_PD);
    alpha_T1=2*tan(alpha_T1);
end

dMTdSPD = @(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT) (ST1.*TRMT.*alpha_PD.*(TRPD.*alpha_T1.^2 - TRT1.*alpha_PD.^2).*(alpha_MT.*SPD.^2.*ST1.*TRPD.^2.*alpha_T1.^2 - alpha_MT.*SPD.^2.*ST1.*TRPD.*TRT1.*alpha_PD.^2 - alpha_MT.*SPD.^2.*ST1.*TRT1.^2.*alpha_T1.^2 + SMT.*SPD.^2.*TRPD.^2.*alpha_T1.*alpha_PD.^2 + 2.*alpha_MT.*SPD.*ST1.^2.*TRT1.^2.*alpha_T1.*alpha_PD - 2.*SMT.*SPD.*ST1.*TRPD.*TRT1.*alpha_T1.^2.*alpha_PD - alpha_MT.*ST1.^3.*TRPD.*TRT1.*alpha_T1.^2 + SMT.*ST1.^2.*TRT1.^2.*alpha_T1.^3))./(2.*SMT.*TRPD.*TRT1.*(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1).^2.*(SPD.*alpha_T1 - ST1.*alpha_PD).^2);

dMTdST1 = @(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT) -(SPD.*TRMT.*alpha_T1.*(TRPD.*alpha_T1.^2 - TRT1.*alpha_PD.^2).*(- alpha_MT.*SPD.^3.*TRPD.*TRT1.*alpha_PD.^2 + 2.*alpha_MT.*SPD.^2.*ST1.*TRPD.^2.*alpha_T1.*alpha_PD + SMT.*SPD.^2.*TRPD.^2.*alpha_PD.^3 - alpha_MT.*SPD.*ST1.^2.*TRPD.^2.*alpha_PD.^2 - alpha_MT.*SPD.*ST1.^2.*TRPD.*TRT1.*alpha_T1.^2 + alpha_MT.*SPD.*ST1.^2.*TRT1.^2.*alpha_PD.^2 - 2.*SMT.*SPD.*ST1.*TRPD.*TRT1.*alpha_T1.*alpha_PD.^2 + SMT.*ST1.^2.*TRT1.^2.*alpha_T1.^2.*alpha_PD))./(2.*SMT.*TRPD.*TRT1.*(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1).^2.*(SPD.*alpha_T1 - ST1.*alpha_PD).^2);

dMTdSMT = @(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT) -(SPD.*ST1.*TRMT.*alpha_MT.*(SPD.*TRT1.*alpha_PD - ST1.*TRPD.*alpha_T1).*(TRPD.*alpha_T1.^2 - TRT1.*alpha_PD.^2))./(2.*SMT.^2.*TRPD.*TRT1.*(SPD.*TRPD.*alpha_PD - ST1.*TRT1.*alpha_T1).*(SPD.*alpha_T1 - ST1.*alpha_PD));

dMT = sqrt( dMTdSPD(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT).^2 .* dSPD.^2 + dMTdST1(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT).^2 .* dST1.^2 + dMTdSMT(SPD,ST1,SMT,alpha_PD,alpha_T1,alpha_MT,TRPD,TRT1,TRMT).^2 .* dSMT.^2);

AdMT     = zeros(size(SPD));
tmp1    = dMT;
tmp1 = max(min(tmp1,threshall.MT),-threshall.MT);
AdMT(dMT>threshall.dMT)     = tmp1(dMT>threshall.dMT);

end
