function hMRI_create_residuals(data,Y,reg,nechoes,NEmap,NSMmap,mpm_params,p,VG,threshall)
% This function creates residuals of the linearized mono-exponential fit 
% per contrast and across contrasts (the latter would be for the ESTATICS model).
% S. Mahammadi 06.09.2019
% modified on 26.08.2020 - "errorESTATICS" map is now in the same units as
% R2s.
% 
% In:
% data          - measured signal
% Y             - modeled signal
% reg           - design matrix
% nechoes       - number of echoes
% NEmap         - nifti
% VG            - target structure of 
% 
% Out:
% maps are written out
% 

    dm = VG.dim;
    Ydata = (reshape(data, [sum(nechoes) prod(dm(1:2))]));
    Ydiff   = exp(Ydata)- exp(reg*Y);
    for ccon = 1:mpm_params.ncon

        Yechotmp = zeros(dm(1:2));
        YMSK = find(sum(Ydata(1+sum(nechoes(1:ccon-1)):sum(nechoes(1:ccon)),:),1)>0);

        if(ccon>1)
            Yechotmp(YMSK) = rms(Ydiff(1+sum(nechoes(1:ccon-1)):sum(nechoes(1:ccon)),YMSK),1);
        else
            Yechotmp(YMSK) = rms(Ydiff(1:nechoes(ccon),YMSK),1);
        end            
        NEmap(ccon).dat(:,:,p) = Yechotmp;
    end
    ccon = ccon + 1;
    tmp = zeros(dm(1:2));
    AR2s = zeros(dm(1:2));
    Ydiff       = Ydata- reg*Y;
    tmp(YMSK)   = rms(Ydiff(:,YMSK),1);
    AR2s(YMSK)        = Y(end,YMSK);
    NEmap(ccon).dat(:,:,p) = tmp*1e3; % in ms
    NSMmap.dat(:,:,p) = (AR2s./tmp).*(tmp>threshall.dR2s); % in ms
    
    
end

