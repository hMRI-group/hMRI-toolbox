function [Yw,w0] = hMRI_FIT_ESTATICS(varargin)
% This function calculates the weighted-least square fit. It requires a
% noise estimate and regularization factor.
% In:
% Y         - ols estimate of data
% DM        - design matrix (reg in MTprot)
% ydata     - log of measured signal
% zpos      - the slice position (p in MTpro)
% thr_w0    - regularization factor
% sigmaMPM  - noise estimate 
%
% Out:
% Yw        - wols fit of data
% w0        - weights
%
% 11/09/2019 S. Mohammadi

for k=1:2:length(varargin),
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% define modelled signal
lSthr   = (DM*Y')';
DDM     = sparse(DM');
% Tikhonov regulization to avoid signularities in large design matrix
LL      = eye(size(DM,2));
LLM     = sparse(LL');
if ~exist('dummyfig','var')
    dummyfig = false;
end

% voxel-wise signal weighting - noise correction for the log-signal
tmp     = ones(size(ydata));
if(sum(sigmaMPM)~=0), tmp     = (bsxfun(@rdivide,exp(lSthr),exp(sigmaMPM))).^2; end

% in case more factors are going to be added to weights 
w0      = tmp;

% normalise to one
mw0 = max(w0,[],2);
MSKtmp = find(mw0>0);
w0(MSKtmp,:) = bsxfun(@rdivide,w0(MSKtmp,:),mw0(MSKtmp));   
% remove all voxels that have smaller weight than zero
minw0 = min(w0,[],2);
maxw0 = max(w0,[],2);
MSKtmp = find(minw0>eps & maxw0<=1);
w0 = w0(MSKtmp,:);   
w0      = w0';
ydata = ydata(MSKtmp,:);
if(dummyfig)
    plot(mean(w0,2),'x-'); ylim([0 1]);title(['Weighting at slice: ' num2str(zpos)]);set(gca,'fontsize',16)

end
www0    = w0(:);
Aw0     = sparse(1:numel(www0),1:numel(www0),www0);
Ew      = speye(size(ydata,1));
X       = kron(Ew,DDM);
yw      = ydata';
yyw     = yw(:);
XX      = (Aw0*X');
XLL       = kron(Ew,LLM);
% regularization
XXX     = X*XX+thr_w0^2*XLL;
% no regularization
%XXX     = X*XX;
clear XX;
YY      = Aw0*yyw;
YYY     = X*YY;
clear YY;
Atmp_w  = XXX\YYY;
Asym_w  = reshape(Atmp_w(:,1),size(DM,2),size(ydata,1));
Asym_w  = Asym_w';

Yw = Y;
Yw(MSKtmp,:) = Asym_w;
Yw(isnan(Asym_w))=Y(isnan(Asym_w));
end
