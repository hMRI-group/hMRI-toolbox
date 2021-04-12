function [parValue, parLocation] = get_metadata_val_classic(varargin)
% This is get_metadata_val, part of the metadata library
%
% PURPOSE
% To retrieve parameter values (mainly acquisition parameters, but this can
% be developed further in the future) from nifti extended header (extended
% nifti format) or simple nifti files with associated JSON metadata. The
% first argument can either be a file name (extended nifti image, normal
% nifti image accompanied by a JSON metadata file, or a JSON metadata file)
% or the matlab structure previously extracted from such files. The second
% argument must be either a predefined parameter name (see list below), or
% a potential match for a given field in the metadata structure. The script
% can therefore be used to efficiently search any (metadata) structure.
%
% IMPORTANT NOTE FOR FUTURE DEVELOPMENTS
% Many cases below are "Siemens-specific", i.e. implemented for Siemens
% data and only valid with Siemens data. Future development will aim at
% making these cases valid for all vendors.
%
% USAGE AND EXAMPLES
% [parValue, parLocation] = get_metadata_val(mstruc, inParName)
% [parValue, parLocation] = get_metadata_val(filenam, inParName)
% parValue = get_metadata_val(filenam, 'RepetitionTime')
%
% ARGUMENTS
% - mstruc is the matlab structure containing the metadata
% - filenam is the name of a file containing or associated with JSON
%   metadata (see above).
% - inParName can be an arbitrary string that will be searched in the
%   metadata structure or one of the following predefined parameter name.
%   Names listed with an '=' are valid for all vendors, while '-' indicate
%   Siemens-only:
%    = 'RepetitionTime' TR [ms]
%    = 'RepetitionTimes' TR [ms]
%    = 'EchoTime' TE(s) [ms]
%    = 'EchoTimes' TE(s) [ms]
%    = 'FlipAngle' [deg]
%    = 'ProtocolName'
%    = 'SequenceName'
%    - 'MT' (1/0 = ON/OFF)
%    = 'FieldStrength' [T]
%    = 'Frequency' [Hz]
%    = 'ScanningSequence'
%    = 'BandwidthPerPixelRO' [Hz/Px]
%    - 'BandwidthPerPixelPE' [Hz/Px]
%    - 'PATparameters' [struct]
%    - 'AccelFactorPE'
%    - 'AccelFactor3D'
%    - 'PELines'
%    - 'PELinesPF'
%    - 'PELinesPAT'
%    - 'PhaseEncodingDirectionPositive' A>>P & R>>L = 1; P>>A & L>>R = 0.
%    = 'PhaseEncodingDirection' A>>P/P>>A = 'COL' & R>>L/L>>R = 'ROW'
%    - 'NumberOfMeasurements'
%    - 'NumberOfSlices'
%    - 'epiReadoutDuration' [ms]
%    - 'WipParameters' structure containing fields alFree & adFree
%    - 'B1mapNominalFAValues' [deg]
%    - 'RFSpoilingPhaseIncrement' [deg]
%    - 'B1mapMixingTime' [ms]
%    - 'AllDiffusionDirections' list of diffusion directions
%    - 'AllBValues' list of b-value
%    - 'DiffusionDirection' diffusion direction of individual DW image
%    - 'BValue' b-value of individual DW image
%    - 'MultiBandFactor' for CMRR multiband and Siemens SMS sequences
%
% - parValue is the value of the parameter. The type of the output can vary
%   according to the request (single real/char/complex value, an array, a
%   cellarray, ...). If inParName is not leading to a unique value, a cell
%   array of values is returned, each element corresponding to a different
%   location in the header structure, as specified by the returned
%   parLocation.
% - parLocation is a string (or a cellarray of strings if the solution is
%   not unique) giving the location of the retrieved value(s) in the
%   metadata structure. For example, if 'RepetitionTime' value is found
%   in the mstruc as mstruc.acqpar.RepetitionTime, parLocation =
%   'acqpar.RepetitionTime'.
%
%==========================================================================
% Written by Evelyne Balteau - June 2016 - Cyclotron Research Centre
%
% December 2016: modified to deal with JSON metadata contained either in
% extended header or separate JSON files.
%==========================================================================

if nargin~=2
    error('Wrong number of arguments. Type: help get_metadata_val for usage.');
end

mstruc = varargin{1};
inParName = varargin{2};
% Should not be accessed directly ?
% % filename as argument, first retrieve the matlab structure:
% if ischar(mstruc)
%     mstruc = get_metadata(mstruc);
%     mstruc = mstruc{1};
% end
% 
parValue = [];
parLocation = [];
nFieldFound = 0;

switch inParName
    case 'RepetitionTime' % [ms]
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'RepetitionTime',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000; % [ms]
        else
          [nFieldFound, fieldList] = find_field_name(mstruc,'RepetitionTimeExitation',...
                                                     'caseSens','sensitive',...
                                                     'matchType','exact');
          [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
          if nFieldFound
              cRes = 1;
              parLocation{cRes} = nam{1};
              parValue{cRes} = val{1} * 1000; % [ms]
          end
        end

        
    case 'RepetitionTimes' % [ms]
        % Not strictly a Bids field
        [nFieldFound, fieldList] = find_field_name(mstruc,'RepetitionTimeSeries',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000; % [ms]
        else
            [parValue, parLocation] = get_metadata_val(mstruc, 'RepetitionTime');
            nFieldFound = size(parLocation, 1);
        end
        
    case 'EchoTime' % [ms]
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'EchoTime',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000; % [ms]
        end
        
    case 'EchoTimes' % [ms]
        % Siemens-specific but made valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'EchoTimeSeries',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000; % [ms]
        else
            [parValue, parLocation] = get_metadata_val(mstruc, 'EchoTime');
            nFieldFound = size(parLocation, 1);
        end
        
    case 'FlipAngle' % [deg]
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'FlipAngle',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    % ProtocolName is not in BIDS json for some reason
        
    case 'SequenceName'
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'SequenceName',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'MT'
        % Not BIDS field
        % may be need to be converted from On/Off
        % if fails get from file name
        [nFieldFound, fieldList] = find_field_name(mstruc,'MTState',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList); %#ok<ASGLU>
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            if ischar(val{1})
              if strcmpi(val{1}, 'on')
                parValue{cRes} = 1;
              elseif strcmpi(val{1}, 'off')
                parValue{cRes} = 0;
              else
                warning(['Invalid MTState value ' val{1}]);
                parValue{cRes} = 0;
              end
            else
              parValue{cRes} = (val{1} > 0);
            end
        else
          % trying from filename
          if isfield(mstruc, 'filename')
            mt = find_entity(mstruc.filename, 'mt-');
            if strcmpi(mt, 'on') || strcmp(mt, '1')
              nFieldFound = 1;
              parLocation{1} = 'filename';
              parValue{1} = 1;
            elseif strcmpi(mt, 'off') || strcmp(mt, '0')
              nFieldFound = 1;
              parLocation{1} = 'filename';
              parValue{1} = 0;
            end
          end
        end
        
    case 'FieldStrength' % [T]
        % Valid for all vendors
        % NB: flNominalB0 returns ~2.8936 for a 3T magnet, but is Siemens
        % specific so unusable with GE/Philips data. However, since more
        % accurate, we first try and retrieve it:
        [nFieldFound, fieldList] = find_field_name(mstruc,'MagneticFieldStrength',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end

    case 'MagneticFieldStrength' % [T]
        % Valid for all vendors
        % NB: flNominalB0 returns ~2.8936 for a 3T magnet, but is Siemens
        % specific so unusable with GE/Philips data. However, since more
        % accurate, we first try and retrieve it:
        [nFieldFound, fieldList] = find_field_name(mstruc,'MagneticFieldStrength',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    % Frequency not BIDS field 
        
    case 'ScanningSequence' % e.g. 'EP' for EPI...
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'ScanningSequence',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
    % BandwidthPerPixelRO not BIDS
    % BandwidthPerPixelPE not BIDS    
    % PELinesPF not BIDS
    % PELines not BIDS
    % PELinesPAT not BIDS
        
    case 'PhaseEncodingDirectionSign'
        % Siemens-specific:
        [nFieldFound, fieldList] = find_field_name(mstruc,'PhaseEncodingDirection',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % PhaseEncodingDirection = [ijk][+-].
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            PEd = val{1};
            if PEd(end) == '-'
              parValue{cRes} = -1;
            else
              parValue{cRes} = 1;
            end
        else
            fprintf(1,['\nWARNING: PhaseEncodingDirection not defined. '...
                     'Using PhaseEncodingDirectionSign directly instead\n']);
            [nFieldFound, fieldList] = find_field_name(mstruc, ...
                            'PhaseEncodingDirectionSign',...
                            'caseSens','sensitive',...
                            'matchType','exact');
          [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
          if nFieldFound
              cRes = 1;
              parLocation{cRes} = nam{1};
              parValue{cRes} = val{1}; 
          end            
        end
        
    case 'PhaseEncodingDirection' % 'COL' (A>>P/P>>A) or 'ROW' (R>>L/L>>R)
        % Valid for all vendors
        [nFieldFound, fieldList] = find_field_name(mstruc,'PhaseEncodingDirection',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            PEd = val{1};
            if PEd(1) == 'j'
              parValue{cRes} = 'COL';
            elseif PEd(1) == 'i'
              parValue{cRes} = 'ROW';
            else
              warning(['Invalid PhaseEncodingDirection ' PEd]);
            end
        end
        
    % NumberOfMeasurements not BIDS
    % NumberOfSlices not BIDS
    % PATparameters not BIDs
    % AccelFactorPE not BIDS (possibly corresponds to MultibandAccelerationFactor)
    % AccelFactor3D not BIDS (possibly corresponds to MultibandAccelerationFactor)
    % MultiBandFactor not BIDS (possibly corresponds to MultibandAccelerationFactor)
    % WipParameters not BIDs

    % AllDiffusionDirections not BIDS (probably equivalent to bvec)
    case 'AllDiffusionDirections'
      if isfield(mstruc, 'filename')
        ext = find_entity(mstruc.filename, 'extension');
        bvec_file = [ mstruc.filename(1:end-size(ext,2)) '.bvec'];
        if exist(bvec_file, 'file')
          nFieldFound = 1;
          parLocation{1} = mstruc.filename;
          parValue{1} = dlmread(bvec_file);
        end
      end

    % AllBValues not BIDS (probably equivalent to bval)    
    case 'AllBValues'
      if isfield(mstruc, 'filename')
        ext = find_entity(mstruc.filename, 'extension');
        bvec_file = [ mstruc.filename(1:end-size(ext,2)) '.bval'];
        if exist(bvec_file, 'file')
          nFieldFound = 1;
          parLocation{1} = mstruc.filename;
          parValue{1} = dlmread(bvec_file);
        end
      end

    % DiffusionDirection not BIDS (can be retrieved from all values with index?)
    % BValue not BIDS (above)

    % isDWI need filename
    case 'isDWI'
      if isfield(mstruc, 'filename')
        nFieldFound = 1;
        suffix = find_entity(mstruc.filename, 'suffix');
        parLocation{1} = 'filename';
        if strcmpi(suffix, 'dwi')
          parValue{1} = 1;
        else
          parValue{1} = 0;
        end

      else
        nFieldFound = 0;
      end
        
    % epiReadout not BIDS but related to TotalReadoutTime
    case 'epiReadoutDuration' % [ms]
        [nFieldFound, fieldList] = find_field_name(mstruc,'epiReadoutTime',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000;
        else
          fprintf(1,['\nWARNING: epiReadoutTime not defined. '...
                     'Using TotalReadoutTime instead\n']);
          [nFieldFound, fieldList] = find_field_name(mstruc,'TotalReadoutTime',...
                                                     'caseSens','sensitive',...
                                                     'matchType','exact');
          [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
          if nFieldFound
              cRes = 1;
              parLocation{cRes} = nam{1};
              parValue{cRes} = val{1} * 1000; 
          end
        end
        
    case 'EchoSpacing' % [ms]
        [nFieldFound, fieldList] = find_field_name(mstruc,'EffectiveEchoSpacing',...
                                                   'caseSens','sensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1} * 1000;
        end
        
    % B1mapNominalFAValues not BIDS, using custom FlipAngleSeries
    case 'B1mapNominalFAValues' % [deg]
        [nFieldFound, fieldList] = find_field_name(mstruc,'FlipAngleSeries',...
                                                   'caseSens','insensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            if isrow(val{1})
              parValue{cRes} = val{1};
            else
              parValue{cRes} = val{1}';
            end
        end
        

    case 'RFSpoilingPhaseIncrement' % [deg]
        [nFieldFound, fieldList] = find_field_name(mstruc,'SpoilingRFPhaseIncrement',...
                                                   'caseSens','insensitive',...
                                                   'matchType','exact');
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
        [nFieldFound, fieldList] = find_field_name(mstruc,'MixingTime',...
                                                   'caseSens','insensitive',...
                                                   'matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
          cRes = 1;
          parLocation{cRes} = nam{1};
          parValue{cRes} = val{1} * 1000; % in ms
        end
end

if ~nFieldFound
    fprintf(1,'\nWARNING: No %s found in the extended header\n', inParName);
    parValue = [];
    parLocation = [];
end

% returns cell array only if necessary (non-unique result)
if (length(parValue) == 1) && iscell(parValue)
    parValue = parValue{1};
    parLocation = parLocation{1};
end
end


function [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList) %#ok<INUSL>
% to retrieve the values and location in the mstruc structure
if nFieldFound
    val = cell(1,nFieldFound);
    nam = cell(1,nFieldFound);
    for cRes = 1:nFieldFound
        nF = 0;
        for cF = 1:size(fieldList,2)
            if ~isempty(fieldList{cRes,cF})
                nF = nF+1;
            end
        end
        nam{cRes} = fieldList{cRes,1};
        for cF = 2:nF
            if iscell(eval(['mstruc.' nam{cRes}]))
                nam{cRes} = [nam{cRes} '{1}'];
                % fprintf(1,['\nWARNING - get_metadata_val/get_val_nam_list:' ...
                %    '\n\tThe value(s) retrieved are one sample out of a bigger cell array.' ...
                %    '\n\tMight be worth checking no other value(s) is(are) available that should be accounted for.\n']);
            elseif length(eval(['mstruc.' nam{cRes}]))>1
                nam{cRes} = [nam{cRes} '(1)'];
                % fprintf(1,['\nWARNING - get_metadata_val/get_val_nam_list:' ...
                %    '\n\tThe value(s) retrieved is(are) one sample out of a bigger array.' ...
                %    '\n\tMight be worth checking no other value(s) is(are) available that should be accounted for.\n']);
            end
            nam{cRes} = [nam{cRes}  '.' fieldList{cRes,cF}];
        end
        val{cRes} = eval(['mstruc.' nam{cRes}]);
    end
else
    val = {};
    nam = {};
end
end

function res = find_entity(fname, field)
  [pth, base, ext] = fileparts(fname);
  base = [base, ext];
  res = [];

  l_field = size(field, 2);
  pos = 0;

  for i = 1:size(base,2)
    if base(i) == '_'
      sub = base(pos + 1: i - 1);
      pos = i;
      if strncmp(sub, field, l_field)
        res = sub(l_field: end);
        break
      end
    elseif base(i) == '.'
      if strcmp('suffix', field)
        res = base(pos + 1: i - 1);
      elseif strcmp('extension', field)
        res = base(i:end);
      end
    end
  end
end
