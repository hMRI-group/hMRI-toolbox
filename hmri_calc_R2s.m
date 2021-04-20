function [R2s,extrapolated]=hmri_calc_R2s(weighted_data)
% R2* estimation using the ESTATICS (ESTimating the Apparent Transverse 
% relaxation time from Images with different ContrastS) model
%
% ledwards@cbs.mpg.de
%
% Based on the ESTATICS implementation in the hMRI toolbox:
%   https://github.com/hMRI-group/hMRI-toolbox
%
% array of structures (one per contrast) input in the form:
%   weighted(contrast).data (NvoxelsX x NvoxelsY x ... x Nechoes)
%   weighted(contrast).TE   (1 x Nechoes)
% Voxels must correspond between the weightings (i.e. the images should 
% have been resliced to the same space), but the sampled TEs may be 
% different.
% Nechoes must be at least 2 for each weighting.
% 
% outputs:
%   R2s (NvoxelsX x NvoxelsY x ...) contains the voxelwise-estimated 
%       common R2* of the weightings.
%   extrapolated contains size(weighted_data) extrapolated to TE=0 in 
%       the same order as the input (e.g. matching contrast order).
% 
% Examples:
%   Ordinary least squares (OLS) estimate of R2* from PD-weighted data:
%       R2s = hmri_calc_R2s(struct('data',PDw,'TE',PDwTE));
%
%   ESTATICS estimate of common R2* of PD, T1, and MT-weighted data, along
%   with weightings extrapolated to TE=0:
%       [R2s,extrapolated] = hmri_calc_R2s([struct('data',PDw,'TE',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TE',MTwTE)]);
%
% Reference:
%     Weiskopf et al. Front. Neurosci., (2014) 
%         https://doi.org/10.3389/fnins.2014.00278

assert(isstruct(weighted_data(1)),'hmri:structError',['inputs must be structs; see help ' mfilename])

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
    
    y=[y;log(rData).']; %#ok<AGROW>
end

%% Solve for R2s
% Well posed as sufficient rank; avoids backslash for data matrix below
DtDSlashDt=(D'*D)\D';

% extra unity in reshape argument to avoid problems if size(dims)==2.
R2s=reshape(DtDSlashDt(1,:)*y,[dims(1:end-1),1]);

%% Extrapolate weightings to TE=0
if nargout>1
    extrapolated=cell(size(weighted_data)); % cell element per contrast
    for w=1:Nweighted
        
        extrapolatedData=exp(DtDSlashDt(w+1,:)*y);
        extrapolated{w}=reshape(extrapolatedData(:),[dims(1:end-1),1]);
        
    end
end

end
