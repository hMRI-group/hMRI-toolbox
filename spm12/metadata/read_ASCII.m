function ascout = read_ASCII(ascin)
%==========================================================================
% function ascout = read_ASCII(ascin)
%
% DESCRIPTION:
% To parse the content of the ASCII part of the DICOM header (Siemens
% specific) and convert it into a matlab structure. Compatible with VA, VB,
% VD and also reading *.SR files (PhoenixProtocols).
%
% USAGE:
% ascout = read_ASCII(ascin)
% where ascin is the char content of either
% hdr.CSASeriesHeaderInfo.MrPhoenixProtocol or
% hdr.CSASeriesHeaderInfo.MrProtocol according to Siemens software version,
% and asc is a matlab structure containing the information contained in the
% ASCCONV BEGIN - ASCCONV END part of the DICOM header. 
%
% NOTE: formatting slightly different from read_ascconv in
% spm_dicom_convert...
%
% WRITTEN BY: Evelyne Balteau - Cyclotron Research Centre and improved
% thanks to the read_ascconv implementation in spm_dicom_convert
%==========================================================================

% only for debugging purpose
ENABLE_DEBUG = false;

% extract the portion between ### ASCCONV BEGIN ### and ### ASCCONV END ###
ascin = regexprep(ascin,'^.*### ASCCONV BEGIN [^#]*###(.*)### ASCCONV END ###.*$','$1');

% split ascin into lines
% asclines = strsplit(ascin,'\n'); % ok for Matlab 2013 and later versions
% replaced by textscan (with delimiter '\n') for compatibility with
% Matlab 2012 and earlier :/...
tmp = textscan(ascin,'%s','delimiter','\n');
asclines = tmp{1};

% do a first cleaning pass over the data:
% - replace indexes [] by () + increment index +1 (C>Matlab indexing)
% - replace double "" and single " by single '
% - replace hexadecimal values (e.g. "0x01") by decimal value
% - delete lines where a field starts or ends by "_"
% - delete end of lines after "#"
asclines = regexprep(asclines,{'\[([0-9]*)\]','["]+','^([^"'']*)0x([0-9a-fA-F]+)',' #.*','^.*\._.*$'},{'($1+1)','''','$1hex2dec(''$2'')','',''});
    
% initialise variables
ascout = [];
clinenum = 1;

while clinenum<length(asclines)
    
    if ENABLE_DEBUG; fprintf(1,'Line #%d - [%s]\n', clinenum, asclines{clinenum}); end %#ok<UNRCH>
    if ~isempty(asclines{clinenum})
        try
            [tlhrh,tmp] = regexp(asclines{clinenum}, '(?:=)+', 'split', 'match'); %#ok<*ASGLU>
            % first process the name of the parameter (hdrnam)
            % split into fields
            [hdrnam,tmp] = regexp(tlhrh{1}, '(?:\.)+', 'split', 'match');
            % check whether any field starts with a number and if so,
            % replace it by "index<number>..."
            isfirstdigit = cellfun(@(x) isstrprop(x(1),'digit'),hdrnam);
            if any(isfirstdigit)
                for cfield=1:length(hdrnam)
                    if isfirstdigit(cfield)
                        hdrnam{cfield} = ['index' hdrnam{cfield}];
                    end
                end
            end
            % concatenate back the fields
            hdrnam = sprintf('.%s', hdrnam{:});
            % now process the value (hdrval)
            hdrval = strtrim(tlhrh{2});
            % and eval the whole line
            eval(sprintf('ascout%s = %s;', hdrnam, hdrval));
            if ENABLE_DEBUG; fprintf(1,'ascout%s = %s;\n', hdrnam, hdrval); end %#ok<UNRCH>
        catch
            fprintf(1,'AscConv: Error evaluating [ %s; ]\n', asclines{clinenum});
        end
    end
    clinenum = clinenum+1;
end
