function tdyhdr = eb_spm_tidycsa(dichdr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% To rearrange CSAImageHeaderInfo and CSASeriesHeaderInfo structures
% obtained when reading the DICOM headers with spm_dicom_headers (SPM12)
% and make it a simplified, easier to browse structure)
%
% USAGE:
% tdyhdr = eb_spm_tidycsa(dichdr)
% where dichdr is a DICOM header structure as read by spm_dicom_headers in
% SPM12.
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
WARNMSGID = 'Matlab:eb_spm_tidycsa';
warning('OFF',WARNMSGID);

if iscell(dichdr);
    dichdr = dichdr{1};
end

tdyhdr = dichdr;

if isfield(dichdr,'CSAImageHeaderInfo')
    tdyhdr = rmfield(tdyhdr,'CSAImageHeaderInfo');
    for ccsa = 1:length(dichdr.CSAImageHeaderInfo)
        val = get_numaris4_val(dichdr.CSAImageHeaderInfo,dichdr.CSAImageHeaderInfo(ccsa).name);
        % if val is empty, let's not waste disk space with empty fields
        if ~isempty(val)
            % if only alphanumeric characters (no letters) we can assume that
            % the value is actually an array of numbers - let's convert them to
            % numbers:
            if (ischar(val) && ~sum(isletter(val(:))))
                tmp = str2num(val);
                if ~isempty(tmp);val = tmp;end
            end
            tdyhdr.CSAImageHeaderInfo.(dichdr.CSAImageHeaderInfo(ccsa).name) = val;
        end
    end
end

if isfield(dichdr,'CSASeriesHeaderInfo')
    tdyhdr = rmfield(tdyhdr,'CSASeriesHeaderInfo');
    for ccsa = 1:length(dichdr.CSASeriesHeaderInfo)
        val = get_numaris4_val(dichdr.CSASeriesHeaderInfo,dichdr.CSASeriesHeaderInfo(ccsa).name);
        % if val is empty, let's not waste disk space with empty fields
        if ~isempty(val)
            % if only alphanumeric characters (no letters) we can assume that
            % the value is actually an array of numbers - let's convert them to
            % numbers:
            if (ischar(val) && ~sum(isletter(val(:))))
                tmp = str2num(val);
                if ~isempty(tmp);val = tmp;end
            end
            tdyhdr.CSASeriesHeaderInfo.(dichdr.CSASeriesHeaderInfo(ccsa).name) = val;
        end
    end
end