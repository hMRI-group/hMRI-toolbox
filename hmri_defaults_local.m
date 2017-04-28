function hmri_defaults_local(varargin)
% PURPOSE
% To set the centre-specific defaults which are used by the hMRI toolbox.
% Customized processing parameters can be defined, overwriting defaults
% from hmri_defaults. Acquisition protocols can also be specified here as a
% fallback solution when no metadata are available. Note that the use of
% metadata is strongly recommended. 
% This script can be conveniently modified by users to fit local settings
% without modifying the global defaults (which is NOT recommended!!). 
%
% WARNING
% Modificaiton of the defaults parameters may impair the the integrity of
% the toolbox, leading to unexpected behaviour. Only recommended for expert
% users. 
%
% HOW DOES IT WORK?
% This local script is called at the end of the standard hmri_defaults
% script. It overwrites and appends the defaults global variable with
% additional fields and values. Additional B1 protocols can be defined and
% made available in the Batch GUI. 
% See examples below.
%_______________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium

%%
global hmri_def

%% ================ definition of all potential protocols =================
if nargin==0 
    
%======================== Global parameters ===============================
% Specifying the lab and scanner
hmri_def.centre = 'crc' ; % 'fil', 'lren', 'crc', 'sciz', 'cbs', ...
hmri_def.scanner = 'prisma' ; % 'allegra', 'terra', ...

%================== B1 mapping processing parameters ======================
% Centre-specific list of protocols available
% NOTE1: the list may contain protocols defined globally and not
%   redefined locally.
% NOTE2: all protocol names MUST 
% - start with a letter, 
% - include only letters, numbers or underscores, i.e. NO space, as these
%   names are used to define a structure fieldname with the protocol
%   parameters. 
% - the first label in the list is the default one.
%==========================================================================
hmri_def.b1_type.labels  = {
    'i3D_EPI'
    'i3D_EPI_local_example'
    'i3D_AFI'
    'tfl_b1_map'
    'rf_map'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    }';
hmri_def.b1_type.val  = hmri_def.b1_type.labels(1);

%==========================================================================
% Now define the above-listed local protocols...
%
% For protocols already defined globally (in hmri_defaults) you may either
% keep the global defaults values, or define your own ones (which will
% overwrite the global ones). In the latter case, proceed as described for
% new protocols below.
%
% For new protocols (see example below), some rules must be followed which
% are described hereafter:
% - The new protocol must be a variation of an existing global one (called
%   template in the following). The variation may be quite extensive.
% - First copy the template parameters.
% - Then modify any field from the template that needs customization.
% - NEVER modify the b1type field!!
% - Good practice for naming a new protocol: <new name> should be the
%   <template name> appended with meaningful suffix (see example below).

%==========================================================================
% 'i3D_EPI' - local example
% In this example, we create a new protocol "i3D_EPI_local_example" which
% is a variation of the i3D_EPI protocol with a few differences redefined
% below. Here we include a modified b1proc parameter as well.
%==========================================================================
hmri_def.b1map.i3D_EPI_local_example = hmri_def.b1map.i3D_EPI;
hmri_def.b1map.i3D_EPI_local_example.b1acq.beta = 140:-7.5:65;
hmri_def.b1map.i3D_EPI_local_example.b1acq.TM = 45;
hmri_def.b1map.i3D_EPI_local_example.b1acq.tert = 19.44; % total EPI RO duration (ms) 
% b0-acquisition
hmri_def.b1map.i3D_EPI_local_example.b0acq.shortTE = 4.72; % ms
hmri_def.b1map.i3D_EPI_local_example.b0acq.longTE = 7.38; % ms
% b1-processing
hmri_def.b1map.i3D_EPI_local_example.b1proc.HZTHRESH = 300;


%% ===================== selecting a specific protocol ====================
% Only defined for selecting default parameters for b1_type, but could be
% extended if needed...
%==========================================================================
elseif nargin==2
    % re-initiate default values
    hmri_defaults;
    % collect arguments
    defstr = varargin{1};
    defval = varargin{2};
    % modify defaults
    switch defstr
        case 'b1_type'
            % retrieve the values
            b1map = hmri_get_defaults(['b1map.' defval]);
            % overwrite global defaults with local ones
            hmri_get_defaults(['b1map.' b1map.b1type],b1map);
    end
    
else 
    error('Wrong number of input! Type "help hmri_defaults_local".');
end
end