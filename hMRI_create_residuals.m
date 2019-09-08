function hMRI_create_residuals(data,Y,reg,nechoes,NEmap,mpm_params,p,VG )
% This function creates residuals of the linearized mono-exponential fit 
% per contrast and across contrasts (the latter would be for the ESTATICS model).
% S. Mahammadi 06.09.2019
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
    tmp(YMSK) = rms(Ydiff(:,YMSK),1);
    NEmap(ccon).dat(:,:,p) = tmp;
    
end

