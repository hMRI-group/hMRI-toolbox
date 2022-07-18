function [R2s,extrapolated,SError]=hmri_calc_R2s(weighted_data,method)
% R2* estimation using an implementation of the ESTATICS
% model (Weiskopf2014). Can utilise weighted least squares (WLS) instead of
% the original ordinary least squares (OLS; Weiskopf2014) to account
% for the heteroscedasticity of log transformed data (Edwards2022).
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
%     -Options are: 'OLS'    (log-linear ordinary least squares estimate),
%                   'WLS[N]' (log-linear weighted least squares estimate
%                             with '[N]' iterations, where '[N]' is 1, 2 ,
%                             or 3; uses OLS signal estimates for initial
%                             weights),
%                   'NLLS_[METHOD]' (non-linear least squares estimate,
%                                    where [METHOD] is the method to be
%                                    used for the initial guess of the
%                                    parameters, e.g. 'ols' or 'wls1'; note
%                                    that while this is expected to be
%                                    accurate, this method is very slow!)
%     -The weights in the WLS case depend on the unknown true signal
%      intensities. These weights can be iteratively updated using the
%      estimated signal intensities. However the benefit of
%      iteratively updating the weights has been found to be relatively
%      small for typical MPM data, and so 'wls1' seems to be sufficient
%      to improve R2* map quality over OLS (Edwards2022).
%
% Outputs:
%   R2s (NvoxelsX x NvoxelsY x ...): the voxelwise-estimated
%       common R2* of the weightings.
%   extrapolated: cell array containing data extrapolated to TE=0 in the
%       same order as the input (e.g. matching contrast order).
%   SError.weighted: Per contrast root mean square error of signal minus 
%       fitted signal. Will be used for error maps.
%   SError.R2s: Root mean square residual of R2* fit.
%
% Examples:
%   OLS estimate of R2* from PD-weighted data array:
%       R2s = hmri_calc_R2s(struct('data',PDw,'TE',PDwTE),'OLS');
%
%   WLS estimate of R2* from PD-weighted data with 3 iterations:
%       R2s = hmri_calc_R2s(struct('data',PDw,'TE',PDwTE),'WLS3');
%
%   NLLS estimate of R2* from PD-weighted data using OLS for initial
%   parameters:
%       R2s = hmri_calc_R2s(struct('data',PDw,'TE',PDwTE),'NLLS_OLS');
%
%   ESTATICS estimate of common R2* of PD, T1, and MT-weighted data,
%   along with estimates of the extrapolation of the input data to TE=0:
%       [R2s,extrapolated] = hmri_calc_R2s([struct('data',PDw,'TEs',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TEs',MTwTE)],'OLS');
%
%   WLS-ESTATICS estimate of common R2* of PD, T1, and MT-weighted data
%   with 1 iteration:
%       R2s = hmri_calc_R2s([struct('data',PDw,'TE',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TE',MTwTE)],'WLS1');
%
% References:
%   Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent
%     transverse relaxation time (R2*) from images with different contrasts
%     (ESTATICS) reduces motion artifacts",
%     https://doi.org/10.3389/fnins.2014.00278
%   Edwards et al. Proc. Int. Soc. Magn. Reson. Med. (2022), "Robust and
%     efficient R2* estimation in human brain using log-linear weighted
%     least squares"
% Mohammadi S, Streubel T, Klock L, Lutti A, Pine K, Weber S, Edwards L,
% Scheibe P, Ziegler G, Gallinat J, Kühn S, Callaghan MF, Weiskopf N,
% Tabelow K (2022)
% Error quantification in multi-parameter mapping facilitates robust
% estimation and enhanced group level sensitivity. :2022.01.11.475846
% Available at: https://www.biorxiv.org/content/10.1101/2022.01.11.475846v1 [Accessed March 5, 2022].

assert(isstruct(weighted_data),'hmri:structError',['inputs must be structs; see help ' mfilename])

dims=size(weighted_data(1).data);
Nvoxels=prod(dims(1:end-1));
Nweighted=numel(weighted_data);

%% Build regression arrays
% Build design matrix
D=[];
for w=1:Nweighted
    d=zeros(length(weighted_data(w).TE),Nweighted+1);
    d(:,1)=-weighted_data(w).TE;
    d(:,w+1)=1;
    D=[D;d]; %#ok<AGROW>
end

% Build response variable vector
y=[];
for w=1:Nweighted
    
    nTEs=length(weighted_data(w).TE);
    assert(nTEs>1,'each weighting must have more than one TE')
    
    localDims=size(weighted_data(w).data);
    assert(localDims(end)==nTEs,'echoes must be in the final dimension')
    assert(prod(localDims(1:end-1))==Nvoxels,'all input data must have the same number of voxels');
    
    rData=reshape(weighted_data(w).data,Nvoxels,nTEs);
    
    % log(0) is not defined, so warn the user about zeroes in their data
    % for methods involving a log transform.
    % The warning can be disabled with "warning('off','hmri:zerosInInput')"
    if any(rData(:)==0)&&~contains(lower(method),'nlls')
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
    
    y=[y;rData.']; %#ok<AGROW>
end

%% Estimate R2*

switch lower(method)
    case 'ols'
        beta=OLS(log(y),D);
        beta(2:end,:)=exp(beta(2:end,:));
    case {'wls1','wls2','wls3'}
        % Number of WLS iterations is specified using 'WLS[N]', where '[N]'
        % is a positive integer
        r=regexp(lower(method),'^wls(\d+)$','tokens');
        niter=str2double(r{1}{1});
        
        logy=log(y);
        
        % Use OLS estimates for initial weights
        y0=exp(D*OLS(logy,D));
        
        % Loop over voxels
        parfor n=1:size(y,2)
            beta(:,n)=WLS(logy(:,n),D,y0(:,n),niter);
        end
        beta(2:end,:)=exp(beta(2:end,:));
        
    case {'nlls_ols','nlls_wls1','nlls_wls2','nlls_wls3'}
        % Check for NLLS case, where specification of the log-linear
        % initialisation method is in the method string following an
        % underscore
        r=regexp(lower(method),'^nlls_(.*)$','tokens');
        initmethod=r{1}{1};
        
        % Initial estimate using log-linear method
        [R2s0,A0]=hmri_calc_R2s(weighted_data,initmethod);
        
        if ~exist('opt','var')
            opt=optimset('lsqcurvefit');
        end
        opt=optimset(opt,'Display','off');
        
        beta0=zeros(Nweighted+1,Nvoxels);
        beta0(1,:)=R2s0(:);
        for w=1:Nweighted
            beta0(w+1,:)=A0{w}(:);
        end
        
        expDecay=@(x,D) (D(:,2:end)*x(2:end)).*exp(x(1)*D(:,1));
        
        % Loop over voxels
        parfor n=1:size(y,2)
            beta(:,n)=lsqcurvefit(@(x,D)expDecay(x,D),beta0(:,n),D,y(:,n),[],[],opt);
        end
    otherwise
        error(['method ' method ' not recognised'])
end

%% Output
% extra unity in reshape argument avoids problems if size(dims)==2.
R2s=reshape(beta(1,:),[dims(1:end-1),1]);

% Extrapolate weightings to TE=0
if nargout>1
    extrapolated=cell(size(weighted_data)); % cell element per contrast
    for w=1:Nweighted
        extrapolatedData=beta(w+1,:);
        extrapolated{w}=reshape(extrapolatedData(:),[dims(1:end-1),1]);
    end
end

% residuals per contrast for error maps
if nargout>2
    
    % output is cell with element per contrast
    SError.weighted=cell(size(weighted_data));
    SError.R2s=zeros([dims(1:end-1),1]);
    for w = 1:Nweighted
        dims=size(weighted_data(w).data);
        TEs=reshape(weighted_data(w).TE,[ones(1,length(dims)-1),dims(end)]);
        
        % per contrast residual
        Ydiff=weighted_data(w).data-extrapolated{w}.*exp(-R2s.*TEs);
        SError.weighted{w} = rms(Ydiff,length(dims)); % rms along TE dimension
    end
    
    % Original implementation of R2* error
    SError.R2s = rms(log(y) - D*[beta(1,:);log(beta(2:end,:))],1);
    SError.R2s = reshape(SError.R2s,[dims(1:end-1),1]);
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
