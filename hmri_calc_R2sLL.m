function [R2s,extrapolated]=hmri_calc_R2sLL(weighted_data,method)
% R2* estimation using an implementation of the ESTATICS
% model (Weiskopf2014) with weighted least squares (WLS) to account for the
% heteroscedasticity of log transformed data (Salvador2005).
%
% ledwards@cbs.mpg.de
%
% Based on the ESTATICS implementation in the hMRI toolbox:
%   https://github.com/hMRI-group/hMRI-toolbox
% and the weighted ESTATICS implementation from Siawoosh Mohammadi:
%   https://github.com/siawoosh/hMRI-toolbox/blob/master/hMRI_FIT_ESTATICS.m
% but unlike Siawoosh's implementation, at present there is no debiasing or
% outlier detection.
% The WLS method applied here differs from the weighted least squares
% implementations used in (Weiskopf2014), which attempted to estimate
% weights based on a method of outlier detection, rather than choosing
% weights to account for the log-transform induced heteroscedasticity.
%
% Input:
%   array of structures (one per contrast) in the form:
%     weighted(contrast).data (NvoxelsX x NvoxelsY x ... x Nechoes)
%     weighted(contrast).TE  (1 x Nechoes)
%   -Voxels must correspond between the weightings (i.e. the images should
%    have been resliced to the same space), but the sampled TEs may be
%    different.
%   -Nechoes must be at least 2 for each weighting.
%   -Because log(0) is ill-defined, zero values in any voxel will result
%    in NaN output for that voxel. To avoid potentially biasing the data,
%    we do not modify the input in any way to avoid this, and leave it to
%    the user to decide how to handle this case, e.g. by removing the
%    corresponding voxels from the input data or replacing zeroes with a
%    small positive number.
%
%   method:
%     string stating which log-linear estimation method to use. 
%     -Options are: 'OLS'    (ordinary least squares estimate), 
%                   'WLS[N]' (weighted least squares estimate with '[N]' 
%                             iterations, where '[N]' is a positive 
%                             integer), and 
%                   'robust' (least squares estimate using matlab's 
%                             'robustfit' function). 
%     -The weights in the WLS case depend on the unknown true signal 
%      intensities. These weights can be iteratively updated using the 
%      estimated signal intensities (Salvador2005). The number of 
%      iterations for convergence is usually 2-3, so this can be a small 
%      number (e.g. use WLS1 or WLS3).
%
% Outputs:
%   R2s (NvoxelsX x NvoxelsY x ...): the voxelwise-estimated
%       common R2* of the weightings.
%   extrapolated: cell array containing data extrapolated to TE=0 in the
%       same order as the input (e.g. matching contrast order).
%
% Examples:
%   OLS estimate of R2* from PD-weighted data:
%       R2s = hmri_calc_R2sLL(struct('data',PDw,'TE',PDwTE),'OLS');

%   WLS estimate of R2* from PD-weighted data with 3 iterations:
%       R2s = hmri_calc_R2sLL(struct('data',PDw,'TE',PDwTE),'WLS3');
%
%   ESTATICS estimate of common R2* of PD, T1, and MT-weighted data,
%   along with estimates of the extrapolation of the input data to TE=0:
%       [R2s,extrapolated] = hmri_calc_R2sLL([struct('data',PDw,'TEs',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TEs',MTwTE)],'OLS');
%
%   WLS-ESTATICS estimate of common R2* of PD, T1, and MT-weighted data
%   with 1 iteration:
%       R2s = hmri_calc_R2sLL([struct('data',PDw,'TE',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TE',MTwTE)],'WLS1');
%
% References:
%   Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent
%     transverse relaxation time (R2*) from images with different contrasts
%     (ESTATICS) reduces motion artifacts",
%     https://doi.org/10.3389/fnins.2014.00278
%   Salvador et al., Hum. Brain Mapp. (2005), "Formal characterization and
%     extension of the linearized diffusion tensor model",
%     https://doi.org/10.1002/hbm.20076

assert(isstruct(weighted_data),'hmri:structError',['inputs must be structs; see help ' mfilename])

dims=size(weighted_data(1).data);
Nvoxels=prod(dims(1:end-1));
Nweighted=numel(weighted_data);

%% Build regression arrays
% Design matrix
D=[];
for w=1:Nweighted
    d=zeros(length(weighted_data(w).TE),Nweighted+1);
    d(:,1)=-weighted_data(w).TE;
    d(:,w+1)=1;
    D=[D;d]; %#ok<AGROW>
end

% Response variable vector
y=[];
for w=1:Nweighted
    
    nTE=length(weighted_data(w).TE);
    assert(nTE>1,'each weighting must have more than one TE')
    
    localDims=size(weighted_data(w).data);
    assert(localDims(end)==nTE,'echoes must be in the final dimension')
    assert(prod(localDims(1:end-1))==Nvoxels,'all input data must have the same number of voxels');
    
    rData=reshape(weighted_data(w).data,Nvoxels,nTE);
    
    % log(0) is not defined, so warn the user about zeroes in their data.
    % The warning can be disabled with "warning('off','hmri:zerosInInput')"
    if any(rData(:)==0)
        warning('hmri:zerosInInput',[...
            'Zero values detected in some voxels in the input data. This ',...
            'will cause estimation to fail in these voxels due to the log ',...
            'transform. If these voxels are background voxels, consider ',...
            'removing them from the input data matrices. ',...
            'Zero values which occur only at high TE in voxels of interest ',...
            'could be replaced with a small positive number, e.g. eps(1) ',...
            '(if the data magnitudes are ~1) or 1 if the data are ',...
            'integer-valued. Note: Care must be taken when replacing ',...
            'values, as this could bias the R2* estimation.']);
    end
    
    y=[y;log(rData).']; %#ok<AGROW>
end

%% Estimate R2*

switch lower(method)
    case 'ols'
        beta=OLS(y,D);
    case 'robust'
        % Loop over voxels
        parfor n=1:size(y,2)
            beta(:,n)=robustfit(D,y(:,n),[],[],'off');
        end
    otherwise
        % Check for WLS case, where the number of WLS iterations is 
        % specified using 'WLS[N]', where '[N]' is a positive integer
        r=regexp(lower(method),'^wls(\d+)$','tokens');
        if ~isempty(r)
            
            niter=str2double(r{1}{1});
            
            % Use OLS estimates for initial weights
            y0=exp(D*OLS(y,D));
            
            % Loop over voxels
            parfor n=1:size(y,2)
                beta(:,n)=WLS(y(:,n),D,y0(:,n),niter);
            end
        else
            error(['method ' method ' not recognised'])
        end
end

%% Output
% extra unity in reshape argument avoids problems if size(dims)==2.
R2s=reshape(beta(1,:),[dims(1:end-1),1]);

% Extrapolate weightings to TE=0
if nargout>1
    extrapolated=cell(size(weighted_data)); % cell element per contrast
    for w=1:Nweighted
        extrapolatedData=exp(beta(w+1,:));
        extrapolated{w}=reshape(extrapolatedData(:),[dims(1:end-1),1]);
    end
end

end

%% Fitting methods
function beta=OLS(y,D)
% allows for vectorized voxel processing

beta=(D'*D)\(D'*y);

end

function beta=WLS(y,D,y0,niter)
% y0 estimate needed for initial weights, could be raw y values or estimate 
% from OLS.
% This function only allows single voxel processing!

% Fix number of iterations to avoid checking convergence of each voxel
for m=1:niter
    % weights are updated using latest parameter estimates
    W=diag(y0.*conj(y0));
    W=W./trace(W); % normalisation
    
    % WLS estimate
    beta=(D'*W*D)\(D'*W*y);
    
    y0=exp(D*beta);
end

end