function quiqi_build=tbx_scfg_hmri_QUIQI_build
% 
% PURPOSE: Compute dictionaries of covariance matrices from the motion 
% degradation index.The dictionary  will be used by spm_reml.m or 
% spm_reml_sc.m to account for heteroscedasticity of the data.
% See " reference to the paper " for more details
% 
%
% METHODS: SPM.mat created when designing the model is modified to add in
% SPM.xVi.Vi the dictionary of covariance matrices. 
%
%_______________________________________________________________________
% Nadege Corbin
% 2021.03.30
% Centre de Resonance Magnetique des Systemes Biologiques, Bordeaux, France
% ======================================================================

% ---------------------------------------------------------------------
% SPM.mat file
% ---------------------------------------------------------------------
spm_mat_file         = cfg_files;
spm_mat_file.tag     = 'spm_mat_file';
spm_mat_file.name    = 'SPM.mat file';
spm_mat_file.help    = {'Select the SPM.mat file containing the design of the model.'};
spm_mat_file.ufilter = '^SPM.mat$';
spm_mat_file.num     = [1 1];

% ---------------------------------------------------------------------
% Power of the MDI included in the dictionary 
% ---------------------------------------------------------------------
lambda        = cfg_entry;
lambda.tag     = 'lambda';
lambda.name    = 'MDI powers ';
lambda.val     = {[0]};
lambda.strtype = 'e';
lambda.num     = [1 Inf];
lambda.help    = {['Specify the powers of the MDI to include in the dictionary. ',... 
    'Vector of integers is expected. By default, a power of 0 is used, '... 
    'equivalent to the OLS case.']};

% ---------------------------------------------------------------------
% MDI json files
% ---------------------------------------------------------------------
MDIfname        = cfg_files;
MDIfname.tag     = 'filename';
MDIfname.name    = 'filename';
MDIfname.help    = {'Select a list of json files containing the quality assessment of the qMRI maps (one file per participant)'};
MDIfname.ufilter = 'hMRI_map_creation_quality_assessment.json$';
MDIfname.num     = [1 Inf];

%----------------------------------------------------------------------
% MDI to be selected from the json file 
% ---------------------------------------------------------------------

MTw        = cfg_menu;
MTw.tag    = 'MTw';
MTw.name   = 'MTw';
MTw.help   = {'MDI from MTw images.'};
MTw.labels = {'yes' 'no'};
MTw.values = {true false};
MTw.val = {false};

T1w        = cfg_menu;
T1w.tag    = 'T1w';
T1w.name   = 'T1w';
T1w.help   = {'MDI from T1w images.'};
T1w.labels = {'yes' 'no'};
T1w.values = {true false};
T1w.val = {false};


PDw        = cfg_menu;
PDw.tag    = 'PDw';
PDw.name   = 'PDw';
PDw.help   = {'MDI from PDw images.'};
PDw.labels = {'yes' 'no'};
PDw.values = {true false};
PDw.val = {true};

%----------------------------------------------------------------------
% Data Type
% ---------------------------------------------------------------------
MDISource = cfg_branch;
MDISource.tag     = 'MDIsource';
MDISource.name    = 'MDI source';
MDISource.help    = {'Use MDI obtained from PD weighted images, T1 weighted ',...
    'images, MT weighted image or a combination.'}; 
MDISource.val = {PDw T1w MTw};


% ---------------------------------------------------------------------
% MDI values 
% ---------------------------------------------------------------------
MDIvalues        = cfg_entry;
MDIvalues.tag     = 'MDIvalues';
MDIvalues.name    = 'MDI values';
MDIvalues.val     = {[ones(10,1) 2*ones(10,1)]};
MDIvalues.strtype = 'e';
MDIvalues.num     = [Inf Inf];
MDIvalues.help    = {['Specify the Motion Degradation Index ',...
    'for each participant (one line per participant). Several columns ',...
    'can be used if several MDI are available (For example, T1 maps have ',...
    '2 MDI available for each participant, one from the PDw images and ',...
    'another one from the T1w images']};


%----------------------------------------------------------------------
% MDI json file
% ---------------------------------------------------------------------

MDIjson          = cfg_branch;
MDIjson.tag       = 'MDIjsonfile';
MDIjson.name      = 'MDI jsonfile';
MDIjson.help      = {'For this type, the user will provide a list of json file.'};
MDIjson.val       = {MDIfname MDISource};

%----------------------------------------------------------------------
% MDI matrix type
% ---------------------------------------------------------------------

MDImatrix           = cfg_branch;
MDImatrix.tag       = 'MDImatrix';
MDImatrix.name      = 'MDI matrix';
MDImatrix.help      = {'For this type, the user will provide a matrix of MDI.'};
MDImatrix.val       = {MDIvalues};


%----------------------------------------------------------------------
% MDI type
% ---------------------------------------------------------------------
MDItype = cfg_choice;
MDItype.tag     = 'MDItype';
MDItype.name    = 'MDI type';
MDItype.help    = {'Type of MDI data. The user can provide a matrix of MDI ',...
    'values or provide a list of json files containing the information.'}; 
MDItype.values  = {MDImatrix MDIjson };
MDItype.val     = {MDImatrix};

% ---------------------------------------------------------------------
% Compute dictionnary of covariance matrices based on the otion degradation
% index
% ---------------------------------------------------------------------

quiqi_build         = cfg_exbranch;
quiqi_build.tag     = 'quiqi_build';
quiqi_build.name    = 'QUIQI BUILD';
quiqi_build.val     = { spm_mat_file lambda MDItype};
quiqi_build.help    = {'Given the MDI index of each participant, a dictionary of '...
    'covariance matrices is built and stored in SPM.xVi.Vi. this dictionary ',...
    'will be used subsequently by spm_reml.m or spm_reml_sc.m to account for ',...
    'the heteroscedasticity of the data when estimating the model parameters.'};
quiqi_build.prog    = @hmri_quiqi_build;
