function hmri_log(varargin)
%==========================================================================
% FORMAT hmri_log(msg,flags)
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
%           formatting. NOTE: the message must be interpreted already, 
%           i.e. must not include \n or \t, hence the recommendation to 
%           use sprintf for the input message!  
%==========================================================================
% PURPOSE
% To be able to review and keep track of info, warnings and other messages
% coming up during processing of the data. The pop-up option is enabled by
% default and can be disabled via the Batch GUI. It is recommended to 
% leave it enabled when piloting the data processing (single subject 
% dataset) to read through and acknowledge every message and make sure
% everything is set up properly before running the processing on a whole
% group.
%
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
    flags.PopUp = false; % message displayed as pop-up window blocking code execution, not enabled by default
end
if ~isfield(flags,'LogFile')
    flags.LogFile = struct('Enabled',true,'LogDir','','FileName',''); % log into file enabled by default
end
if ~isfield(flags.LogFile,'Enabled')
    error('Flags not properly set to call hmri_log, check format');
end

% uninterpreted \t and \n must be converted before escaping the \ character
% msg = strrep(msg,'\t',sprintf('\t'));
% msg = strrep(msg,'\n',sprintf('\n'));
% NO: this might be problematic when printing a directory path with leading
% "t" and "n" in directory names :/... The interpretation must be done when
% calling hmri_log: always use "hmri_log(sprintf(...), flags)"...

% for fprintf, some special characters must be escaped:
fpfmsg = msg;
fpfmsg = strrep(fpfmsg,'\','\\');
fpfmsg = strrep(fpfmsg,'%','%%');

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
    fprintf(fid,['\n' fpfmsg '\n' ]);
    fclose(fid);
end

if flags.ComWin
    fprintf(1,['\n' fpfmsg '\n' ]);
end

if flags.PopUp
    % creade message window in modal mode (to restrict user interaction to
    % that specific window):
    h = msgbox(msg,'hMRI toolbox','modal');
    
    ah = get( h, 'CurrentAxes' );
    ch = get( ah, 'Children' );
    
    % change font size and font name for readability
    fs = 12;
    fn = 'Courier';
    set(ch, 'FontName', fn, 'FontSize', fs);
    
    % readjust window size
    textExtent = get(ch, 'Extent');
    winPos = get(h,'Position');
    textOffset = 7;
    OKbuttonOffset = 20;
    set(h, 'Position', [winPos(1) winPos(2) textExtent(3)+textOffset*2 textExtent(4)+textOffset*2+OKbuttonOffset]); 

    % use uiwait to block the code execution until msgbox acknowledged
    uiwait(h);
end


end
       

