%=========================================================================%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%=========================================================================%
% Copyright (C) 2014 - Cyclotron Research Centre
% University of Liege, Belgium

function hdr = eb_read_ASCII(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% this script retrieves protocol parameters from the ASCII section in the
% DICOM header. Version compatible with VD13D and reading from PHOENIX *.SR
% files.
%
% USAGE:
% hdr = eb_read_ASCII(filename)
% where filename is the header file name (format meas*.hdr) and hdr is a
% structure containing all the retrieved parameters.
%
% WARNING AND DISCLAIMER:
% This software is for research use only. Do not use it for clinical or
% diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - November 2014
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENABLE_DEBUG = false;

fid = fopen(filename);
if fid == -1
    disp('File does not exist or cannot be opened');
    return;
end

% initialise a few variables
hdr = [];
% cdepth = 1;
new_field_name = '';

% starts reading...
%tic
cline = fgets(fid); %returns the next line of the file, if end-of-file, returns -1
clinenum = 1;
isasciihdr = 0;

numbers = '1234567890';

while (cline ~= -1)
    if ENABLE_DEBUG
        disp(['Line #' num2str(clinenum) ' - [' cline(1:end-1) ']']);
    end
    
    idx_asconv_begin = strfind(cline,'### ASCCONV BEGIN'); % for Prisma data
    idx_asconv_end = strfind(cline,'### ASCCONV END ###');
    
    if ~isempty(idx_asconv_begin)
        disp('Reading ASCII headers');
        isasciihdr = 1;
        new_field_name = 'asc';
        if isfield(hdr, new_field_name)
            rep = eval(['length(hdr.' new_field_name ')'])+1;
        else
            rep = 1;
        end
        new_field_name = [new_field_name '{' num2str(rep) '}'];
        cline = fgets(fid); % go to next line where ascii headers start
    end
    if ~isempty(idx_asconv_end)
        isasciihdr = 0;
        
        %break;
    end
    
    if isasciihdr
        [hdrnam, hdrval] = strtok(cline,'=');
        hdrval = strtok(hdrval,'=');
        hdrval = hdrval(1:end-1); % to remove the end-of-line character...
        hdrval = strtrim(hdrval); % to remove leading and trailing white space from string
        
        % first process hdrnam
        % skip if contains '_' characters (usually parameter attribute
        % in Prisma data)
        if isempty(strfind(hdrnam,'_'))
            % convert indexes if any: C++ [0],[1],... -> matlab (1),(2),...
            idx = [strfind(hdrnam,'[');strfind(hdrnam,']')];
            if ~isempty(idx)
                % remplace [] by ()
                hdrnam(idx(1,:)) = '('; hdrnam(idx(2,:)) = ')';
                % increment indexes
                for in = 1:(size(idx,2))
                    oldidx = hdrnam(idx(1,in)+1:idx(2,in)-1);
                    newidx = num2str(str2double(oldidx)+1);
                    hdrnam = [hdrnam(1:idx(1,in)) newidx hdrnam(idx(2,in):end)];
                    if (length(newidx)>length(oldidx) && in<size(idx,2))
                        idx(:,in+1:end) = idx(:,in+1:end)+1;
                    end
                end
            end
            % check whether any field starts with a number
            idx = strfind(hdrnam, '.');
            okidx = idx;
            for cidx = 1:length(idx)
                if strfind(numbers,hdrnam(idx(cidx)+1))
                    okidx(cidx) = 0;
                else
                    okidx(cidx) = 1;
                end
            end
            if ~isempty(find(okidx==0,1))
                idx = idx(okidx==0);
                for cidx = length(idx):-1:1
                    hdrnam = [hdrnam(1:idx(cidx)) 'index' hdrnam(idx(cidx)+1:end)];
                end
            end
            % now process hdrval
            idx = strfind(hdrval, '#'); % remove trailing comments (following #)
            if ~isempty(idx)
                hdrval = strtrim(hdrval(1:idx(cidx)-1));
            end
            idx = strfind(hdrval, '""'); % strings are surrounded by ""
            if ~isempty(idx)
                if ENABLE_DEBUG 
                    disp(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval(idx(1)+2:idx(end)-1) '''};']);
                end
                eval(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval(idx(1)+2:idx(end)-1) '''};']);
            elseif ~isempty(strfind(hdrval, 'x')) 
                % hexadecimal numbers stored as strings '0x...'
                if ENABLE_DEBUG
                    disp(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval '''};']);
                end
                eval(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval '''};']);
%                 % hexadecimal numbers are converted into decimal (the '0x'
%                 % prefix must be removed for conversion)
%                 idx = strfind(hdrval, 'x');
%                 if ENABLE_DEBUG 
%                     disp(['hdr.' new_field_name '.' hdrnam ' = hex2dec(''' hdrval(idx+1:end) ''');']);
%                 end
%                 eval(['hdr.' new_field_name '.' hdrnam ' = hex2dec(''' hdrval(idx+1:end) ''');']);
            else
                if ENABLE_DEBUG 
                    disp(['hdr.' new_field_name '.' hdrnam ' = str2num(hdrval);']);
                end
                eval(['hdr.' new_field_name '.' hdrnam ' = str2num(hdrval);']);
            end
        end
    end
    cline = fgets(fid);
    clinenum = clinenum+1;
    
end
if isfield(hdr,'asc')
    if length(hdr.asc)==1
        hdr = hdr.asc{1};
    else
        hdr = hdr.asc{1};
        fprintf('\nWARNING: several ASCII header found in file \n\t%s',filename);
    end
else
    hdr = [];
    fprintf('\nWARNING: no ASCII header found in file \n\t%s',filename);
end
fclose(fid);
%toc