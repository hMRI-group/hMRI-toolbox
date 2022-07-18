function dPD = hmri_make_dPD(SPD,ST1,dSPD,dST1,alpha_PD,alpha_T1,TRPD,TRT1,f_T,small_angle_approximation)
% Calculate propagation of uncertainty for PD map
% (https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
% % S.Mohammadi 06.09.2019
%
% In:
% SPD           - PDw signal at TE=0
% ST1           - T1w signal at TE=0
% dSPD          - residual of mono-exponential fit of PDw signal
% dST1          - residual of mono-exponential fit of T1w signal
% alpha_PD      - flip angle of PDw signal
% alpha_T1w     - flip angle of T1w signal
% TRPD          - repetition time of PDw signal
% TRT1          - repetition time of T1w signal
% f_T           - map of transmit field
%
% Out:
% dPD           - error for A in [a.u.]

if(~isempty(f_T))
    alpha_PD = alpha_PD.*f_T;
    alpha_T1 = alpha_T1.*f_T;
    if ~small_angle_approximation
        alpha_PD=2*tan(alpha_PD/2);
        alpha_T1=2*tan(alpha_T1/2);
    end
end

% dPD calculation is symmetric with respect to the two weighted contrasts
dPD = sqrt( dPD_by_dS1(SPD,ST1,alpha_PD,alpha_T1,TRPD,TRT1).^2.*dSPD.^2 ...
    +dPD_by_dS1(ST1,SPD,alpha_T1,alpha_PD,TRT1,TRPD).^2.*dST1.^2 );

end

function d = dPD_by_dS1(S1,S2,alpha1,alpha2,TR1,TR2)
% Derivative of dual flip-angle A (PD) estimate with respect to first 
% weighted signal (S1). Because of symmetry in the R1 calculation, the 
% derivative with respect to the second weighted signal can be computed by 
% permuting labels.
%
% Can be derived using: 
%   syms S1 alpha1 S2 alpha2
%   syms TR1 TR2 positive
%   diff(hmri_calc_A(struct('data',S1,'fa',alpha1,'TR',TR1,'B1',1),struct('data',S2,'fa',alpha2,'TR',TR2,'B1',1),true),S1)

d = S1.*S2.*TR2.*alpha1.*(TR1.*alpha2./alpha1 - TR2.*alpha1./alpha2) ./ (S1.*TR2.*alpha1 - S2.*TR1.*alpha2).^2 ...
    - S2.*(TR1.*alpha2./alpha1 - TR2.*alpha1./alpha2)./(S1.*TR2.*alpha1 - S2.*TR1.*alpha2);

end
