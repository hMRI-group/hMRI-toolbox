function dicom = spm_cfg_dicom
% SPM Configuration file for DICOM Import
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_dicom.m 6952 2016-11-25 16:03:13Z guillaume $

% ---------------------------------------------------------------------
% data DICOM files
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'DICOM files';
data.help    = {'Select the DICOM files to convert.'};
data.filter  = 'any';
data.ufilter = '.*';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% root Directory structure for converted files
% ---------------------------------------------------------------------
root         = cfg_menu;
root.tag     = 'root';
root.name    = 'Directory structure for converted files';
root.help    = {'Choose root directory of converted file tree. The options are:'
                ''
                '* Output directory: ./<StudyDate-StudyTime>/<ProtocolName>'
                '* Output directory: ./<PatientID>/<ProtocolName>'
                '* Output directory: ./<PatientID>/<StudyDate-StudyTime>/<ProtocolName>'
                '* Output directory: ./<ProtocolName>'
                '* No directory hierarchy: Convert all files into the output directory, without sequence/series subdirectories'}';
root.labels  = {'Output directory: ./<StudyDate-StudyTime>/<ProtocolName>'
                'Output directory: ./<PatientID>/<ProtocolName>'
                'Output directory: ./<PatientID>/<StudyDate-StudyTime>/<ProtocolName>'
                'Output directory: ./<ProtocolName>'
                'No directory hierarchy'}';
% removed 'Output directory: ./<PatientName>/<ProtocolName>' for anonymity purposes
root.values  = {'date_time'
                'patid'
                'patid_date'
                'series'
                'flat'}';
root.def     = @(val)spm_get_defaults('dicom.root', val{:});
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where files are written.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% protfilter Protocol name filter
% ---------------------------------------------------------------------
protfilter         = cfg_entry;
protfilter.tag     = 'protfilter';
protfilter.name    = 'Protocol name filter';
protfilter.help    = {'A regular expression to filter protocol names. DICOM images whose protocol names do not match this filter will not be converted.'};
protfilter.strtype = 's';
protfilter.num     = [0 Inf];
protfilter.val     = {'.*'};
% ---------------------------------------------------------------------
% format Output image format
% ---------------------------------------------------------------------
format         = cfg_menu;
format.tag     = 'format';
format.name    = 'Output image format';
format.help    = {'DICOM conversion can create separate img and hdr files or combine them in one file. The single file option will help you save space on your hard disk, but may be incompatible with programs that are not NIfTI-aware.'
                  'In any case, only 3D image files will be produced.'}';
format.labels  = {'Two file (img+hdr) NIfTI'
                  'Single file (nii) NIfTI'}';
format.values  = {'img' 'nii'};
format.def     = @(val)spm_get_defaults('images.format', val{:});
% ---------------------------------------------------------------------
% JSON metadata format
% ---------------------------------------------------------------------
mformat         = cfg_menu;
mformat.tag     = 'mformat';
mformat.name    = 'JSON metadata format';
mformat.help    = {'Metadata (acquisition parameters, processing steps, ...) ' ...
                    'can be stored using JSON format. The JSON metadata can ' ...
                    'either (or both) be stored as a separate JSON file ' ...
                    'or as an extended header (included in the image file, ' ...
                    'for nii output format only).'}';
mformat.labels  = {'None'
                    'Separate JSON file'
                    'Extended nii header'
                    'Both'}';
mformat.values  = {'' 'sep' 'ext' 'sepext'};
mformat.val     = {''};

% ---------------------------------------------------------------------
% Options for metadata content: anonymisation
% NOTE: anonymisation CANNOT BE GUARANTEED!!
% ---------------------------------------------------------------------
manonym         = cfg_menu;
manonym.tag     = 'manonym';
manonym.name    = 'Anonymisation (attempt!)';
manonym.help    = {['If JSON metadata are to be stored, you might want to ' ...
                   'make sure that patient confidentiality is preserved. '] ...
                   ['WARNING: effective anonymisation cannot be guaranteed, ' ...
                   'since it depends on the way patient data have been ' ...
                   'registered. USE WITH CARE! '] ...
                   'Two levels of anonymisation are available:' ...
                   [' - full anonymisation: no patient information is kept at all.' ...
                   ' NOTE: the patient ID is kept in the converted file name.'] ...
                   [' - basic anonymisation: patient ID (assuming it does not ' ...
                   'contain his name), age (years at the time of the data ' ...
                   'acquisition), sex, size and weight are kept, patient ' ...
                   'name, date of birth and DICOM filename (often containing ' ...
                   'the patient name) are removed.']}';
manonym.labels  = {'None'
                   'Full'
                   'Basic'}';
manonym.values  = {'none','full','basic'};
manonym.val     = {'basic'};
% ---------------------------------------------------------------------
% Options for metadata content: full or essentials
% ---------------------------------------------------------------------
mcontent         = cfg_menu;
mcontent.tag     = 'mcontent';
mcontent.name    = 'Content of JSON metadata';
mcontent.help    = {'To reduce storage space (but is it still a concern?) ' ...
                       'metadata can be limited to a basic set of parameters ' ...
                       'or not...'}';
mcontent.labels  = {'Full'
                    'Essentials'}';
mcontent.values  = {1 0};
mcontent.val     = {1};
% ---------------------------------------------------------------------
% Options for metadata content
% ---------------------------------------------------------------------
metaopts       = cfg_branch;
metaopts.tag   = 'metaopts';
metaopts.name  = 'JSON metadata';
% NOTE: for release version, the anonymisation and content options are
% removed from the Batch...
% metaopts.val   = {mformat mcontent manonym};
metaopts.val   = {mformat};
metaopts.help  = {''};
% ---------------------------------------------------------------------
% icedims Use ICEDims in filename
% ---------------------------------------------------------------------
icedims        = cfg_menu;
icedims.tag    = 'icedims';
icedims.name   = 'Use ICEDims in filename';
icedims.help   = {'If image sorting fails, one can try using the additional SIEMENS ICEDims information to create unique filenames. Use this only if there would be multiple volumes with exactly the same file names.'};
icedims.labels = {'No' 'Yes'};
icedims.values = {0 1};
icedims.val    = {0};
% ---------------------------------------------------------------------
% convopts Conversion options
% ---------------------------------------------------------------------
convopts       = cfg_branch;
convopts.tag   = 'convopts';
convopts.name  = 'Conversion options';
convopts.val   = {format icedims metaopts};
convopts.help  = {''};
% ---------------------------------------------------------------------
% dicom DICOM Import
% ---------------------------------------------------------------------
dicom          = cfg_exbranch;
dicom.tag      = 'dicom';
dicom.name     = 'DICOM Import';
dicom.val      = {data root outdir protfilter convopts};
dicom.help     = {
    'DICOM Conversion.'
    'Most scanners produce data in DICOM format. This routine attempts to convert DICOM files into SPM compatible image volumes, which are written into the current directory by default. Note that not all flavours of DICOM can be handled, as DICOM is a very complicated format, and some scanner manufacturers use their own fields, which are not in the official documentation at http://medical.nema.org/'
    }';
dicom.prog     = @spm_run_dicom;
dicom.vout     = @vout;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job) %#ok<INUSD>
dep            = cfg_dep;
dep.sname      = 'Converted Images';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
