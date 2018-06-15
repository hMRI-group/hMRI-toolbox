function hmri_log(varargin)
%==========================================================================
% FORMAT hmri_log(flags,msg)
%   flags   ComWin: messages logged to the Matlab Command Window;
%           PopUp: messages appearing in a pop-up window, interrupting
%               the processing until acknowledge by the user - useful
%               for piloting data processing and checking every step;
%           LogFile: 
%               Enabled: messages printed in a log file
%               LogDir: directory where to save the log file
%                       (default: current directory (pwd))
%               FileName: log file name (default: hMRI_LogFile.log)
%   msg     Message to be logged, typically created with sprintf for
%           formatting etc...
%==========================================================================
% PURPOSE
% To be able to review and keep track of info, warnings nad other messages
% coming up during processing of the data. The pop-up option is enabled in
% the defaults and can be disabled via the Configure toolbox module. It is
% recommended to leave it enabled when piloting the data processing (single
% subject) to read through and acknowledge every message and make sure
% everything is set up properly before running the processing on a whole
% group.
%
% TEMPORARY NOTE (June 2018): the pop-up option is still at the development
% level. Will be included in the toolbox soon...
%==========================================================================
% Written by E. Balteau, 2018.
% Cyclotron Research Centre, University of Liege, Belgium
%==========================================================================
msg = '';
flags = struct;
if nargin==0
    error('Not enough input arguments');
end
if nargin>=1
    msg = varargin{1};
end
if nargin==2
    flags = varargin{2};
end

if ~isfield(flags,'ComWin')
    flags.ComWin = true; % log into command window enabled by default
end
if ~isfield(flags,'PopUp')
    flags.PopUp = false; % currently in development, not enabled
end
if ~isfield(flags,'LogFile')
    flags.LogFile = struct('Enabled',true,'LogDir','','FileName',''); % log into file enabled by default
end
if ~isfield(flags.LogFile,'Enabled')
    error('Flags not properly set to call hmri_log, check format');
end

if flags.LogFile.Enabled
    if ~isfield(flags.LogFile,'LogDir')
        flags.LogFile.LogDir = '';
    end
    if ~isfield(flags.LogFile,'FileName')
        flags.LogFile.FileName = '';
    end
    if isempty(flags.LogFile.LogDir)
        flags.LogFile.LogDir = pwd;
    end
    if isempty(flags.LogFile.FileName)
        flags.LogFile.FileName = 'hMRI_LogFile.log';
    end
    logfnam = fullfile(flags.LogFile.LogDir,flags.LogFile.FileName);
    fid = fopen(logfnam,'a');
    fprintf(fid,['\n' msg '\n' ]);
    fclose(fid);
end

if flags.ComWin
    fprintf(1,['\n' msg '\n' ]);
end

if flags.PopUp
    % creade message window in modal mode (to restrict user interaction to
    % that specific window):
    h = msgbox(msg,'hMRI toolbox','modal');
    % use uiwait to block the code execution until msgbox acknowledged
    uiwait(h);
end


end
       

