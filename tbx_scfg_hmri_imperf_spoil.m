function imperf_spoil=tbx_scfg_hmri_imperf_spoil
% 
% PURPOSE: Compute correction factors for imperfect spoiling based on the 
% method described in (Preibisch & Deichmann, MRM, 2009). 
%
% METHODS: Numerical simulations are performed with the Extended Phase 
% Graph framework(code from https://github.com/mriphysics/EPG-X). Sequence 
% parameters, a range of values for B1+ efficiency and T1, and global 
% parameters T2 and diffusion coefficent T must be provided. 
%
%_______________________________________________________________________
% Wellcome Centre for Human Neuroimaging
% Nadège Corbin - May 2020
% ======================================================================

% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where output json file will be written to.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% Name of the protocol
% ---------------------------------------------------------------------
prot_name        = cfg_entry;
prot_name.tag     = 'prot_name';
prot_name.name    = 'Protocol Name ';
prot_name.val     = {'UnitTest_protocol'};
prot_name.strtype = 's';
prot_name.help    = {['Specify the name of the protocol']};
% ---------------------------------------------------------------------
% T1 range [ms]
% ---------------------------------------------------------------------
T1range        = cfg_entry;
T1range.tag     = 'T1range_ms';
T1range.name    = 'T1 range ';
T1range.val     = {[500:100:2000]};
T1range.strtype = 'e';
T1range.num     = [1 Inf];
T1range.help    = {['Specify the range of T1 values over which the ',...
    'corrections factors will be computed. A linear fitting ',...
    'T1=A*T1app+B will be performed to estimate A and B for each B1+ '...
    'efficiency.']};

% ---------------------------------------------------------------------
% B1+ efficiency range [%]
% ---------------------------------------------------------------------
B1range        = cfg_entry;
B1range.tag     = 'B1range';
B1range.name    = 'B1 range ';
B1range.val     = {[0.7 : 0.05 : 1.3]};
B1range.strtype = 'e';
B1range.num     = [1 Inf];
B1range.help    = {['Specify the range of transmit field efficiency (B1) over which the ',...
    'corrections factors will be computed. After the linear fitting ',...
    'T1=A*T1app+B , a polynomial fitting will be perfomed to estimate ',...
    'A and B such as: A=P(B1) and B=P(B1) ',...
    ' with P , 2nd degree polynom. '...
    'Note: B1=1 corresponds to optimal transmit field efficiency']};

% ---------------------------------------------------------------------
% T2 [ms] 
% ---------------------------------------------------------------------
T2        = cfg_entry;
T2.tag     = 'T2_ms';
T2.name    = 'T2';
T2.val     = {[70]};
T2.strtype = 'r';
T2.num     = [1 1];
T2.help    = {['Specify an estimate of the T2 in ms']};


% ---------------------------------------------------------------------
% D [um^2/ms] 
% ---------------------------------------------------------------------
D        = cfg_entry;
D.tag     = 'D_um2_per_ms';
D.name    = 'D';
D.val     = {[0.8]};
D.strtype = 'r';
D.num     = [1 1];
D.help    = {['Specify an estimate of the diffusion coeffcient (D) in um^2/ms']};

% ---------------------------------------------------------------------
% Readout gradient amplitude [ms] 
% ---------------------------------------------------------------------
Gdur        = cfg_entry;
Gdur.tag     = 'Gdur_ms';
Gdur.name    = 'Readout gradient duration';
Gdur.val     = {[1.9998]};
Gdur.strtype = 'e';
Gdur.num     = [1 Inf];
Gdur.help    = {['Specify the duration (in ms) of the gradient in the readout direction',...
                'of the FLASH acquisitions ']};
% ---------------------------------------------------------------------
% Readout gradient amplitude [ms] 
% ---------------------------------------------------------------------
Gamp        = cfg_entry;
Gamp.tag     = 'Gamp_mT_per_m';
Gamp.name    = 'Readout gradient amplitude';
Gamp.val     = {[44.04]};
Gamp.strtype = 'e';
Gamp.num     = [1 Inf];
Gamp.help    = {['Specify the amplitude (in mT/m) of the gradient in the readout direction',...
                'of the FLASH acquisitions ']};

% ---------------------------------------------------------------------
% TR [ms] 
% ---------------------------------------------------------------------
TR         = cfg_entry;
TR.tag     = 'TR_ms';
TR.name    = 'TR';
TR.val     = {[25]};
TR.strtype = 'e';
TR.num     = [1];
TR.help    = {['Specify the TR (in ms) of the FLASH acquisitions']};


% ---------------------------------------------------------------------
% RF Spoiling increment [deg]
% ---------------------------------------------------------------------
Phi0         = cfg_entry;
Phi0.tag     = 'Phi0_deg';
Phi0.name    = 'RF spoiling increment';
Phi0.val     = {[137]};
Phi0.strtype = 'w';
Phi0.num     = [1 1];
Phi0.help    = {['Specify the RF SPoiling increment (in deg) of the FLASH acquisitions']};

% ---------------------------------------------------------------------
% Flip angles [deg]
% ---------------------------------------------------------------------
FA        = cfg_entry;
FA.tag     = 'FA_deg';
FA.name    = 'Flip angles';
FA.val     = {[6 21]};
FA.strtype = 'e';
FA.num     = [1 2];
FA.help    = {['Specify the flip angles (in deg) of the PD-weighted ',...
    'and the T1-weighted FLASH acquisitions (in that order)']};

% ---------------------------------------------------------------------
% All tissue parameters
% ---------------------------------------------------------------------
tissue_params            = cfg_branch;
tissue_params.tag        = 'tissue_params';
tissue_params.name       = 'Tissue parameters';
tissue_params.help       = {'Input all the tissue parameters.'};
tissue_params.val        = {T1range T2 D};

% ---------------------------------------------------------------------
% All sequence parameters
% ---------------------------------------------------------------------
seq_params            = cfg_branch;
seq_params.tag        = 'seq_params';
seq_params.name       = 'Sequence parameters';
seq_params.help       = {'Input all the sequence parameters.'};
seq_params.val        = {FA TR Phi0 B1range Gdur Gamp};

% ---------------------------------------------------------------------
% Compute correction factors for imperfect spoiling
% ---------------------------------------------------------------------

imperf_spoil         = cfg_exbranch;
imperf_spoil.tag     = 'imperf_spoil';
imperf_spoil.name    = 'Correction for imperfect spoiling';
imperf_spoil.val     = { outdir prot_name seq_params tissue_params };
imperf_spoil.help    = {'hMRI map creation based on multi-echo FLASH sequences including optional receive/transmit bias correction.'};
imperf_spoil.prog    = @hmri_corr_imperf_spoil;
