function hmri_quiqi_build(job)
%==========================================================================
%
% PURPOSE: Compute dictionaries of covariance matrices from the motion
% degradation index.The dictionary  will be used by spm_reml.m or
% spm_reml_sc.m to account for heteroscedasticity of the data.
% See " reference to the paper " for more details
%
%
% METHODS: SPM.mat created when designing the model is modified to add in
% SPM.xVi.Vi the dictionary of covariance matrices.
%_______________________________________________________________________
% Antoine Lutti
% 2021.04.01
% Neuroimaging Research Laboratory, Lausanne University Hospital &
% University of Lausanne, Lausanne, Switzerland
% Nadege Corbin
% 2021.03.30
% Centre de Resonance Magnetique des Systemes Biologiques, Bordeaux, France
%==========================================================================

hmri_log(sprintf('\t--- Build dictionary of covariance matrices based on the MDI ---'));
%% ***********************************************%%
% Get Inputs
%*************************************************%%

% SPM.mat file
spm_mat_file      =   job.spm_mat_file;
load(spm_mat_file{1});

% get MDI values
MDItype        =   job.MDItype;

if isfield(MDItype,'MDIjsonfile')
    ListFile=MDItype.MDIjsonfile.filename;
    
    if length(ListFile)~= size (SPM.xX.X,1)
        error('The number of json file is different from the number of lines in the design matrix')
    end
    
    
    MDIsource= MDItype.MDIjsonfile.MDIsource;
    for f=1:length(ListFile)
        QAdata=spm_jsonread(ListFile{f});
        col=0;
        if (MDIsource.PDw)
            if isfield(QAdata.SDR2s,'PDw')
                col=col+1;
                MDIvalues(f,col)=QAdata.SDR2s.PDw;
            else
                error('There is no MDI value from PDw data')
            end
        end
        if (MDIsource.T1w)
            if isfield(QAdata.SDR2s,'T1w')
                col=col+1;
                MDIvalues(f,col)=QAdata.SDR2s.T1w;
            else
                error('There is no MDI value from T1w data')
            end
        end
        if (MDIsource.MTw)
            if isfield(QAdata.SDR2s,'MTw')
                col=col+1;
                MDIvalues(f,col)=QAdata.SDR2s.MTw;
            else
                error('There is no MDI value from MTw data')
            end
        end
    end
    
    
else
    MDIvalues      =   job.MDItype.MDImatrix.MDIvalues;
    
    if size(MDIvalues,1)~= size (SPM.xX.X,1)
        error('The number of lines in the MDI matrix is different from the number of lines in the design matrix')
    end
end


% vector of powers of MDI
lambda = job.lambda;

%% ***********************************************%%
% create vectors of MDI to the power of lambda
% ************************************************%%
MatCovDict=[];
for indMDI= 1:size(MDIvalues,2)
    for indLam=1:size(lambda,2)
        MatCovDict=cat(2,MatCovDict,MDIvalues(:,indMDI).^lambda(indLam));
    end
end


%% ***********************************************%%
% Check if a dictionary is already present.
% One will exist if a group comparison has been specified
%*************************************************%%

if isfield(SPM.xVi,'Vi')
    % check if it is a group comparison
    nVi=size(SPM.xVi.Vi,2);
    for v=1:nVi
        len(v)=length(find(diag(SPM.xVi.Vi{v})~=0));
    end
    if length(unique(len))>1  % This is a group comparison with different sizes for both groups 
        [Uni inda indc]=unique(len);
        GroupIndx={};
        for ctr=1:length(unique(len))
            GroupIndx{ctr}=find(diag(SPM.xVi.Vi{inda(ctr)})~=0);
        end
    else
        if unique(len)==size(SPM.xVi.Vi{1},1)/2 % This is a comparison of two groups with the same size 
            GroupIndx{1}=find(diag(SPM.xVi.Vi{1})~=0);
            GroupIndx{2}=find(diag(SPM.xVi.Vi{end})~=0);
        else % This is not a group comparison 
            GroupIndx{1}=find(diag(SPM.xVi.V)~=0);
        end
    end
else % this is not a group comparison 
    GroupIndx{1}=find(diag(SPM.xVi.V)~=0);
end
SPM=rmfield(SPM,'xVi');


%% ***********************************************%%
% Create dictionary of covariance matrices and
% store it in SPM.xVi.Vi
%*************************************************%%
ind=0;
for indGroup=1:size(GroupIndx,2)% separate basis function for each power of the MDI AND group
    for indMat=1:size(MatCovDict,2)
        ind=ind+1;
        DiagTerms=zeros(length(MDIvalues),1);
        DiagTerms(GroupIndx{indGroup},1)=MatCovDict(GroupIndx{indGroup},indMat);
        SPM.xVi.Vi(ind)={sparse(diag(DiagTerms))};
    end
end


%% ***********************************************%%
% Add MDI vaalues to the SPM structure
%*************************************************%%
SPM.QUIQI_MDI=MDIvalues; 

%% ***********************************************%%
% Save the changes in SPM.mat
%*************************************************%%
[pn fn]= fileparts(spm_mat_file{1});
save(fullfile(pn,'SPM.mat'), 'SPM')

end
