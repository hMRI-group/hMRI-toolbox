function varargout = vbq_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defaults = vbq_get_defaults
% Return the global "defaults" variable defined in vbq_defaults.m.
%
% FORMAT defval = vbq_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr".
% Currently, this is a '.' subscript reference into the global
% "vbq_def" variable defined in vbq_defaults.m.
%
% FORMAT vbq_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of VBQ. To make
% persistent changes, edit vbq_defaults.m.
%
% NOTE, specific to VBQ tool:
% In order to allow centre specific defaults, *without* editing/commenting
% the default file itself, these centre specific values are placed in a
% centre specific substructure, named with 'fil', 'lren' or 'crc'(see the 
% vbq_defaults.m file).
% Then if the required field is not found in the standard default 
% structure, the centre specific field, e.g. 'fil', 'lren' or crc', is 
% included *automatically* in the call. Therefore you can access a centre
% specific parameter by specifying the centre properly in the defaults
% file AND calling for the parameter.
% 
% Example:
% --------
% Definition of vbq_def
% vbq_def.centre = 'crc' ; 
% vbq_def.param1 = 123 ;
% vbq_def.crc.TR  = 3;  % in sec
% vbq_def.fil.TR  = 2;  % in sec
% vbq_def.lren.TR  = 2.5;  % in sec
% 
% v = vbq_get_defaults('param1')
% returns the value 123 into v
% v = vbq_get_defaults('TR')
% returns the value 3 into v
% If you edit the file and set vbq_def.centre = 'fil' ;
% then v = vbq_get_defaults('TR')
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
% Then modified for use with the VBQ toolbox by Christophe Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: vbq_get_defaults.m 37 2013-12-18 13:30:38Z christophe $

global vbq_def;
if isempty(vbq_def)
    vbq_defaults;
end

if nargin == 0
    varargout{1} = vbq_def;
    return
end

try
    % Assume it's working as standard SPM functionality
    % construct subscript reference struct from dot delimited tag string
    tags = textscan(defstr,'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    
    if nargin == 1
        varargout{1} = subsref(vbq_def, subs);
    else
        vbq_def = subsasgn(vbq_def, subs, varargin{1});
    end
catch %#ok<CTCH>
    % Try adding the centre name as intermediate field
    ctr_name = vbq_get_defaults('centre');
    % construct subscript reference struct from dot delimited tag string
    tags = textscan([ctr_name,'.',defstr],'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    if nargin == 1
        varargout{1} = subsref(vbq_def, subs);
    else
        vbq_def = subsasgn(vbq_def, subs, varargin{1});
    end
end

end

%% Some demo stuff
% %
% % get the defaults in a standard routine:
% % - a single parameter
% v = vbq_get_defaults('param1')
% v2 = vbq_get_defaults('set1.prefix')
% 
% % - a set of parameters (substructure)
% s = vbq_get_defaults('set1')
% 
% % - for the default centre, one parameter
% vc = vbq_get_defaults('TR')
% 
% % - for the default centre, one set of parameters
% sc = vbq_get_defaults('cset2')
% %
% % in the batch system, use the following syntax
% % - a centre specific parameter from the 'cset1' set
% name.param    = @(val)vbq_get_defaults('cset1.param1', val{:});

