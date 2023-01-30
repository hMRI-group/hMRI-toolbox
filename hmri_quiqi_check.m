function hmri_quiqi_check(job)
%==========================================================================

% PURPOSE: Plot Residuals of the estimated model with respect to the Motion
% Degradation Index to evaluate the efficiency of the weighting
%
%
% METHODS: Spatial variance of the residuals of each subject is plotted
% with respect to their corresponfing MDI value
%
%_______________________________________________________________________
% Antoine Lutti
% 2021.04.01
% Neuroimaging Research Laboratory, Lausanne University Hospital &
% University of Lausanne, Lausanne, Switzerland
% Nadege Corbin
% 2021.03.30
% Centre de Resonance Magnetique des Systemes Biologiques, Bordeaux, France
%==========================================================================

hmri_log(sprintf('\t--- Evaluate residuals and MDI relationship ---'));
%% ***********************************************%%
% Get Inputs
%*************************************************%%

% SPM.mat file
spm_mat_file      =   job.spm_mat_file;
load(spm_mat_file{1});

% Order of the polynomial fit
pow      =   job.power;


%% ***********************************************%%
% Get Residuals files
%*************************************************%%
[pn fn]=fileparts(spm_mat_file{1});
ResFiles=cellstr(spm_select('FPList',pn,'^Res_'));

if length(ResFiles{1})==0
    error('No residual files could be found in the folder of the SPM.mat file. Please ensure to save the image residuals when estimating the model. ')
end


%% ***********************************************%%
% Get Mask
%*************************************************%%

fileMask=cellstr(spm_select('FPList',pn,'^mask.nii'));

if length(fileMask{1})==0
    error('No mask could be found in the folder of the SPM.mat file')
end

Mask=spm_read_vols(spm_vol(fileMask{1}));
MaskIndx=find(Mask~=0);

%% ***********************************************%%
% Get MDI values
%*************************************************%%

if ~(isfield(SPM,'QUIQI_MDI'))
    error('There is no QUIQI_MDI field in the SPM structure')
end

MDIVals=SPM.QUIQI_MDI;

%% ***********************************************%%
% Compute variance of the residuals
%*************************************************%%

ResidVar=zeros(size(MDIVals,1),1);
for Subjctr=1:size(MDIVals,1)% reads-in individual residual maps and estimates variance
    
    tempRes=spm_read_vols(spm_vol(ResFiles{Subjctr}));
    ResidVar(Subjctr)=var(tempRes(MaskIndx),'omitnan');
    
end


%% ***********************************************%%
% Fit of the residuals with respect to MDI
%*************************************************%%
% ResidVar=ResidVar*1e6;% for R2s maps in ms-1 only
[~,Rsq,yfit]=myPolyFit(MDIVals,ResidVar,pow,'Free');
plotResFit(ResidVar,yfit,Rsq,40,SPM.swd)
    
end

function [P,Rsq,yfit]=myPolyFit(X,Y,Powers,FitMethod)
% Polynomial fitting routine
% INPUTS:
%     - X: independent variable
%     - Y: dependent variable
%     - Powers: order of the polynomial fit
%     - FitMethod: allowed values: 1. 'NonNeg' - enforces positive
%     polynomial coefficients. 2. 'Free' - allows positive and negative
%     polynomial coefficients.

% OUTPUTS:
%     - P: polynomial coefficients
%     - Rsq: R-square of the fit
%     - yfit: fitted dependent variable
%
Powers=linspace(0,Powers,Powers+1);
Xmat=[];
for ctr2=1:size(Powers,2)
    if Powers(ctr2)==0
        Xmat=cat(2,Xmat,ones(size(X,1),1));
    else
        Xmat=cat(2,Xmat,X.^Powers(ctr2));
    end
end
if strcmp(FitMethod,'NonNeg')
    P = lsqnonneg(Xmat,Y);
elseif strcmp(FitMethod,'Free')
    P=pinv(Xmat'*Xmat)*Xmat'*Y;
end
yfit =Xmat * P;
SSresid = sum((Y - yfit).^2);
SStotal = (length(Y)-1) * var(Y);
Rsq = 1 - SSresid/SStotal;

end


function plotResFit(Y,yfit,Rsq,Nbins,SavePath)
% Plots residuals against the fit of the residuals with a polynominal
% function of the MDI

figure
hold on
plot(Y,yfit,'*')
plot(Y,Y,'r','Linewidth',2)%for reference
title(['R^2 = ' num2str(round(Rsq*1e2)/1e2)])
ylabel('Residual fit');xlabel('Residuals');
saveas(gcf, fullfile(SavePath,'ResFitvsRes'), 'fig');

end
