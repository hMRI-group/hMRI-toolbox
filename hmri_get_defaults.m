function varargout = hmri_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defaults = hmri_get_defaults
% Return the global "defaults" variable defined in hmri_defaults.m.
%
% FORMAT defval = hmri_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr".
% Currently, this is a '.' subscript reference into the global
% "hmri_def" variable defined in hmri_defaults.m.
%
% FORMAT hmri_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of hMRI. To make
% persistent changes, edit hmri_defaults.m.
%
% NOTE, specific to hMRI tool (might become obsolete in future versions):
% In order to allow centre specific defaults, *without* editing/commenting
% the default file itself, these centre specific values are placed in a
% centre specific substructure, named with 'fil', 'lren' or 'crc'(see the 
% hmri_defaults.m file).
% Then if the required field is not found in the standard default 
% structure, the centre specific field, e.g. 'fil', 'lren' or crc', is 
% included *automatically* in the call. Therefore you can access a centre
% specific parameter by specifying the centre properly in the defaults
% file AND calling for the parameter.
% 
% Example:
% --------
% Definition of hmri_def
% hmri_def.centre = 'crc' ; 
% hmri_def.param1 = 123 ;
% hmri_def.crc.TR  = 3;  % in sec
% hmri_def.fil.TR  = 2;  % in sec
% hmri_def.lren.TR  = 2.5;  % in sec
% 
% v = hmri_get_defaults('param1')
% returns the value 123 into v
% v = hmri_get_defaults('TR')
% returns the value 3 into v
% If you edit the file and set hmri_def.centre = 'fil' ;
% then v = hmri_get_defaults('TR')
% returns the value 2 into v
%
% A few constraints for this trick to work:
% - the default structure has a field called 'centre' defining the current
%   centre values to use ('fil', crc', 'lren' to begin with)
% - all centre substructures should have the same organisation!
% - centre specific fields and global fields should have *different* names,
%   e.g. do NOT define "def.TE = 20" and "def.fil.TE = 30" fields as the 
%   latter would NEVER be used.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% Then modified for use with the hMRI toolbox by Christophe Phillips
% Cyclotron Research Centre, University of Liege, Belgium

global hmri_def;
if isempty(hmri_def)
    hmri_defaults;
end

if nargin == 0
    varargout{1} = hmri_def;
    return
end

try
    % Assume it's working as standard SPM functionality
    % construct subscript reference struct from dot delimited tag string
    tags = textscan(defstr,'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    
    if nargin == 1
        varargout{1} = subsref(hmri_def, subs);
    else
        hmri_def = subsasgn(hmri_def, subs, varargin{1});
    end
catch %#ok<CTCH>
    % Try adding the centre name as intermediate field
    ctr_name = hmri_get_defaults('centre');
    % construct subscript reference struct from dot delimited tag string
    tags = textscan([ctr_name,'.',defstr],'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    if nargin == 1
        varargout{1} = subsref(hmri_def, subs);
    else
        hmri_def = subsasgn(hmri_def, subs, varargin{1});
    end
end

end

%% Some demo stuff
% %
% % get the defaults in a standard routine:
% % - a single parameter
% v = hmri_get_defaults('param1')
% v2 = hmri_get_defaults('set1.prefix')
% 
% % - a set of parameters (substructure)
% s = hmri_get_defaults('set1')
% 
% % - for the default centre, one parameter
% vc = hmri_get_defaults('TR')
% 
% % - for the default centre, one set of parameters
% sc = hmri_get_defaults('cset2')
% %
% % in the batch system, use the following syntax
% % - a centre specific parameter from the 'cset1' set
% name.param    = @(val)hmri_get_defaults('cset1.param1', val{:});

