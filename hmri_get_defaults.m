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
% persistent changes, see help section in hmri_defaults.m and
% hmri_local_defaults.m.
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
    varargout{1} = [];
    fprintf(1,'WARNING: no default value defined for %s!\n', defstr);
end

end