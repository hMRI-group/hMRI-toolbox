function asc = eb_read_phoenix(phoenix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% To parse the content of hdr.CSASeriesHeaderInfo.MrPhoenixProtocol (which
% contains the ASCII part of the DICOM header) and convert it into an
% easier to handle matlab structure.
%
% USAGE:
% asc = eb_read_phoenix(phoenix)
% where phoenix is the char content of
% hdr.CSASeriesHeaderInfo.MrPhoenixProtocol, and hdr is the header read by
% spm_dicom_headers and tidied by eb_spm_tidycsa.  
% asc is a matlab structure containing the information contained in the
% ASCCONV BEGIN - ASCCONV END part of the DICOM header. 
%
% WARNING AND DISCLAIMER:
% This software is for research use only. Do not use it for clinical or
% diagnostic purposes.
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
% Written by Evelyne Balteau - June 2015
% Copyright (C) 2015 - Cyclotron Research Centre
% University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global WARNMSGID ENABLE_DEBUG;
ENABLE_DEBUG = false;
WARNMSGID = 'Matlab:eb_read_phoenix';
warning('OFF',WARNMSGID);

% In order to use eb_read_ASCII code, the txt content of MrPhoenixProtocol
% must be written into a temporary file (will be deleted afterwards).
fid = fopen('tmpASCCONV.txt','w');
fprintf(fid,'%s',phoenix);
fclose(fid);

% Now read the file
asc = eb_read_ASCII('tmpASCCONV.txt');

% And tidy up
delete('tmpASCCONV.txt');