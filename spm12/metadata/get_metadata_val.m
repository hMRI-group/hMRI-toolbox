function [parValue, parLocation] = get_metadata_val(varargin)
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
%    - 'epiReadoutDuration' [ms]
%    - 'WipParameters' structure containing fields alFree & adFree
%    - 'B1mapNominalFAValues' [deg]
%    - 'RFSpoilingPhaseIncrement' [Hz]
%    - 'B1mapMixingTime' [ms]
%    - 'AllDiffusionDirections' list of diffusion directions
%    - 'AllBValues' list of b-value
%    - 'DiffusionDirection' diffusion direction of individual DW image
%    - 'BValue' b-value of individual DW image
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

% filename as argument, first retrieve the matlab structure:
if ischar(mstruc)
    mstruc = get_metadata(mstruc);
    mstruc = mstruc{1};
end

parValue = [];
parLocation = [];
nFieldFound = 0;

if ~isstruct(mstruc)
    % there seems to be no metadata available (mstruc is not a structure
    % and is likely to be empty). Still try a few tricks in case either
    % TR/TE/FA is the target field and available in the description field:
    if ischar(varargin{1}) && strcmp(spm_file(varargin{1},'ext'),'nii')
        fnam = varargin{1};
        fprintf(1,'WARNING: No metadata available in %s.\n', fnam);
        N = nifti(fnam);
        p = struct('tr',[],'te',[],'fa',[]);
        tmp = regexp(N.descrip,'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg','names');
        if ~isempty(tmp)
            p.tr = str2double(tmp.tr);
            p.te = str2double(tmp.te);
            p.fa = str2double(tmp.fa);
            switch inParName
                case 'RepetitionTime' % [ms]
                    cRes = 1;
                    parLocation{cRes} = 'NiftiDescriptionField';
                    parValue{cRes} = p.tr;
                    fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f ms\n', inParName, parValue{1});
                    
                case 'EchoTime'
                    cRes = 1;
                    parLocation{cRes} = 'NiftiDescriptionField';
                    parValue{cRes} = p.te;
                    fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f ms\n', inParName, parValue{1});
                    
                case 'FlipAngle'
                    cRes = 1;
                    parLocation{cRes} = 'NiftiDescriptionField';
                    parValue{cRes} = p.fa;
                    fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f ms\n', inParName, parValue{1});
                    
                otherwise
                    fprintf(1,'Not able to retrieve any value for parameter %s.\n', inParName);
            end
        end
    else
        error('Invalid input arguments. Type help get_metadata_val for correct syntax.');
    end
    
else
    switch inParName
        case 'RepetitionTime' % [ms]
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'RepetitionTime', 'caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'RepetitionTimes' % [ms]
            % Siemens-specific but made valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'alTR', 'caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % alTR is a Siemens-specific field which holds an array of
            % values rather than a single value -> convenient when several
            % TRs used for a given sequence (e.g. AFI). If not available
            % (GE or Philips), let's get the standard DICOM field
            % RepetitionTime instead (in ms):
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1}*0.001;
            else
                [parValue, parLocation] = get_metadata_val(mstruc, 'RepetitionTime');
            end
            
        case 'EchoTime' % [ms]
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'EchoTime', 'caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'EchoTimes' % [ms]
            % Siemens-specific but made valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'alTE', 'caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % alTE is a Siemens-specific field which holds an array of TE
            % values rather than a single value -> convenient if several
            % TEs used for a given sequence (e.g. multiecho sequences). If
            % not available (GE or Philips), let's get the standard DICOM
            % field EchoTime instead (in ms):
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1}*0.001;
            else
                [parValue, parLocation] = get_metadata_val(mstruc, 'EchoTime');
            end
            
        case 'FlipAngle' % [deg]
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'FlipAngle', 'caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many - no scaling necessary
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'ProtocolName'
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'ProtocolName','caseSens','sensitive','matchType','exact');
            if nFieldFound==0 % may happen when Anonymous data, no exact match for ProtocolName (tProtocolName instead)
                [nFieldFound, fieldList] = find_field_name(mstruc, 'ProtocolName','caseSens','sensitive');
            end
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % to correct for mismatch of outputs depending on data
            % Anonymous or not ('tProtocolName' returns a cell array of
            % cell array of char, while 'ProtocolName' returns a cell array
            % of char): 
            if iscell(val) && (nFieldFound==1)
                if iscell(val{1})
                    val = val{1};
                end
            end
            
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'SequenceName'
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'SequenceName','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'MT'
            % Siemens-specific
            % NB: parameters set to 0 are usually omitted in the DICOM header.
            % Therefore, the absence of the parameter in the header means that
            % no MT pulse is applied. If applied, parameter
            % acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPrepPulses.ucMTC
            % takes the hexadecimal value '0x1' = 1.
            [nFieldFound, fieldList] = find_field_name(mstruc, 'ucMTC','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList); %#ok<ASGLU>
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = 1;
            else
                nFieldFound = 1;
                cRes = 1;
                parLocation{cRes} = 'NullParameterNotDefinedInStruct';
                parValue{cRes} = 0;
            end
            
        case 'FieldStrength' % [T]
            % Valid for all vendors
            % NB: flNominalB0 returns ~2.8936 for a 3T magnet, but is Siemens
            % specific so unusable with GE/Philips data. However, since more
            % accurate, we first try and retrieve it:
            [nFieldFound, fieldList] = find_field_name(mstruc, 'flNominalB0','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            else
                % NB: MagneticFieldStrength returns 3 for a 3T magnet
                [nFieldFound, fieldList] = find_field_name(mstruc, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                if nFieldFound
                    cRes = 1;
                    parLocation{cRes} = nam{1};
                    parValue{cRes} = val{1};
                end
            end
            
            
        case 'Frequency' % [Hz]
            % Valid for all vendors
            % NB: lFrequency returns 123255074 Hz and is Siemens-specific,
            % while ImagingFrequency returns 123.2551 MHz.
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lFrequency','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            else
                [nFieldFound, fieldList] = find_field_name(mstruc, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                if nFieldFound
                    cRes = 1;
                    parLocation{cRes} = nam{1};
                    parValue{cRes} = val{1}*1000000;
                end
            end
            
        case 'ScanningSequence' % e.g. 'EP' for EPI...
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'ScanningSequence','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'BandwidthPerPixelRO'
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'PixelBandwidth','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'BandwidthPerPixelPE' % Siemens-specific header entry
            [nFieldFound, fieldList] = find_field_name(mstruc, 'BandwidthPerPixelPhaseEncode','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'PELinesPF'
            % size of the k-space PE dimension, taking into account partial
            % Fourier but not Parallel acceleration. Valid for all vendors.
            [nFieldFound, fieldList] = find_field_name(mstruc, 'NumberOfPhaseEncodingSteps','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'PELines'
            % size of the k-space PE dimension, without taking into account
            % Parallel acceleration nor partial Fourier.
            % Siemens-specific towards validation for all vendors...
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lPhaseEncodingLines','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            else
                % This variation is an attempt for compatibility for all
                % vendors, but not tested thoroughly.
                % PEdir = get_metadata_val(mstruc, 'InPlanePhaseEncodingDirection');
                PEdir = get_metadata_val(mstruc, 'PhaseEncodingDirection');
                if strcmp(deblank(PEdir),'ROW');
                    fieldName = 'Rows';
                else
                    fieldName = 'Columns';
                end
                [parValue, parLocation] = get_metadata_val(mstruc,fieldName);
                if ~isempty(parValue);nFieldFound = 1;end
            end
            
        case 'PELinesPAT'
            % size of the k-space PE dimension, taking into account Parallel
            % acceleration but not partial Fourier. Used to calculate the total
            % EPI Readout duration for FieldMap undistortion.
            % Siemens-specific.
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lPhaseEncodingLines','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            [nFieldFoundPAT, fieldListPAT] = find_field_name(mstruc, 'lAccelFactPE','caseSens','sensitive','matchType','exact');
            [valPAT,namPAT] = get_val_nam_list(mstruc, nFieldFoundPAT, fieldListPAT); %#ok<ASGLU>
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                if isempty(valPAT);valPAT{1} = 1;end
                parValue{cRes} = floor(val{1}/valPAT{1});
            end
            
        case 'PhaseEncodingDirectionSign'
            % Siemens-specific:
            [nFieldFound, fieldList] = find_field_name(mstruc, 'PhaseEncodingDirectionPositive','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % PhaseEncodingDirectionPositive = 0/1 for -1/+1.
            % Note that null parameters are not saved in the header.
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            else
                nFieldFound = 1;
                cRes = 1;
                parLocation{cRes} = 'NullParameterNotDefinedInStruct';
                parValue{cRes} = -1;
            end
            
        case 'PhaseEncodingDirection' % 'COL' (A>>P/P>>A) or 'ROW' (R>>L/L>>R)
            % Valid for all vendors
            [nFieldFound, fieldList] = find_field_name(mstruc, 'InPlanePhaseEncodingDirection','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'NumberOfMeasurements' % Siemens-specific
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lRepetitions','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1}+1;
            else
                nFieldFound = 1;
                cRes = 1;
                parLocation{cRes} = 'NullParameterNotDefinedInStruct';
                parValue{cRes} = 1;
            end
            
        case 'PATparameters'
            % Siemens-specific
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sPat','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'AccelFactorPE'
            % Siemens-specific
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lAccelFactPE','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'AccelFactor3D'
            % Siemens-specific
            [nFieldFound, fieldList] = find_field_name(mstruc, 'lAccelFact3D','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'WipParameters'
            % Siemens-specific (NB: search made case insensitive since
            % sWiPMemBlock or sWipMemBlock depending on software version)
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','insensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
            
        case 'AllDiffusionDirections'
            % Siemens-specific
            if get_metadata_val(mstruc,'isDWI')
                [nFieldFound, fieldList] = find_field_name(mstruc, 'sDiffusion','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                % sDiffusion is the field containing series diffusion information.
                % Example:
                % mstruc{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sDiffusion
                %           lDiffWeightings: 2
                %            lNoiseLevel: 30
                %        lDiffDirections: 117
                %                 ulMode: 128
                %               dsScheme: 2
                %       ulQSpaceCoverage: 1
                %       ulQSpaceSampling: 1
                %       lQSpaceMaxBValue: 50
                %           lQSpaceSteps: 1
                %               alBValue: [0 2000]
                %             alAverages: [1 1]
                %     sFreeDiffusionData: [1x1 struct]
                % with
                % mstruc{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sDiffusion.sFreeDiffusionData
                %               sComment: [1x1 struct]
                %        lDiffDirections: 117
                %     ulCoordinateSystem: 1
                %        ulNormalization: 1
                %        asDiffDirVector: {1x116 cell}
                if nFieldFound
                    % check it is a DWI sequence
                    if isfield(val{1},'sFreeDiffusionData')
                        cRes = 1;
                        % NB: the length of asDiffDirVector may be shorter than the
                        % number of directions since 0 array values are discarded
                        % from the DICOM header (here, the last "direction" is a b0
                        % scan with direction (0,0,0)).
                        ndir = eval(['mstruc.' nam{1} '.sFreeDiffusionData.lDiffDirections']);
                        parLocation{cRes} = [nam{1} '.sFreeDiffusionData.asDiffDirVector'];
                        parValueSagCorTra = eval(['mstruc.' nam{1} '.sFreeDiffusionData.asDiffDirVector']);
                        parValue{cRes} = zeros(3,ndir);
                        for cdir = 1:length(parValueSagCorTra)
                            if isempty(parValueSagCorTra(cdir).dSag);parValueSagCorTra(cdir).dSag = 0;end
                            if isempty(parValueSagCorTra(cdir).dCor);parValueSagCorTra(cdir).dCor = 0;end
                            if isempty(parValueSagCorTra(cdir).dTra);parValueSagCorTra(cdir).dTra = 0;end
                            parValue{cRes}(:,cdir) = [parValueSagCorTra(cdir).dSag; parValueSagCorTra(cdir).dCor; parValueSagCorTra(cdir).dTra];
                        end
                    end
                    
                end
            end
            
        case 'AllBValues'
            % Siemens-specific
            if get_metadata_val(mstruc,'isDWI')
                [nFieldFound, fieldList] = find_field_name(mstruc, 'alBValue','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                if nFieldFound
                    cRes = 1;
                    parLocation{cRes} = nam{1};
                    parValue{cRes} = val{1};
                end
            end
            
        case 'DiffusionDirection'
            % Siemens-specific
            if get_metadata_val(mstruc,'isDWI')
                [nFieldFound, fieldList] = find_field_name(mstruc, 'DiffusionGradientDirection','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                if nFieldFound
                    cRes = 1;
                    parLocation{cRes} = nam{1};
                    parValue{cRes} = val{1};
                else
                    % B0 images have "direction" set to [0 0 0].
                    % Null parameters in DICOM headers are often omitted,
                    % therefore no "direction" information is retrievable for
                    % these images. If non-existent field
                    % "DiffusionGradientDirection" found in DWI acquisition,
                    % vector [0 0 0] is now returned.
                    nFieldFound = 1;
                    cRes = 1;
                    parLocation{cRes} = 'B0Image';
                    parValue{cRes} = [0;0;0];
                    fprintf(1,'WARNING: Diffusion direction not defined for DWImage %s. Assuming b=0.\n', inParName);
                end
            end
            
        case 'BValue'
            % Siemens-specific
            if get_metadata_val(mstruc,'isDWI')
                [nFieldFound, fieldList] = find_field_name(mstruc, 'B_value','caseSens','sensitive','matchType','exact');
                [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
                if nFieldFound
                    cRes = 1;
                    parLocation{cRes} = nam{1};
                    parValue{cRes} = val{1};
                end
            end
            
        case 'isDWI'
            % Siemens-specific
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sDiffusion','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = [nam{1} '.ulMode'];
                if (val{1}.ulMode>1)
                    parValue{cRes} = 1;
                else
                    parValue{cRes} = 0;
                end
            else
                nFieldFound = 1;
                cRes = 1;
                parLocation{cRes} = 'NotDWISequence';
                parValue{cRes} = 0;
            end
            if parValue{1}==0
                fprintf(1,'WARNING: This is not a DWI sequence.\n');
            end
            
            
        case 'epiReadoutDuration' % [ms]
            % Siemens-specific
            % This information is easily retrievable from standard EPI
            % sequences, where the "BandwidthPerPixelPhaseEncoding" is defined.
            % For the 3D-EPI B1mapping sequence, everything is hard coded in
            % the sequence and not passed to the MrProt variable, therefore it
            % is not available in the DICOM header. We have to work case by
            % case, relying on sequence version, BWPP and resolution :/...
            
            % first check whether BandwidthPerPixelPhaseEncode is defined
            [nFieldFound, fieldList] = find_field_name(mstruc, 'BandwidthPerPixelPhaseEncode','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                parValue{cRes} = 1/val{1}*1000;
            else
                fprintf(1,['WARNING: BandwidthPerPixelPhaseEncode not defined for the current sequence\n' ...
                    'For 3D-EPI B1 mapping sequences, values are deduced from the sequence\n' ...
                    'version. Be aware that it might not be correct if the version is unknown.\n']);
                
                % check whether it is an EPI sequence (al_B1mapping is not defined as 'EP' but 'RM'):
                valEPI = get_metadata_val(mstruc, 'ScanningSequence');
                valSEQ = get_metadata_val(mstruc, 'SequenceName');
                valPROT = get_metadata_val(mstruc, 'ProtocolName');
                
                % Case al_B1mapping
                if strfind(lower(valPROT),'b1map')
                    fprintf(1,'Trying to derive the epiReadoutDuration from the Sequence version (%s/%s).\n',valSEQ,valPROT);
                    nFieldFound = 1;
                    switch lower(valSEQ)
                        case 'b1v2d3d2' % 540 us - Prisma
                            EchoSpacing = 2*140+260; 
                        case 'b1epi4a3d2' % 330 us - Allegra
                            EchoSpacing = 330; 
                        case 'b1epi2b3d2' % 1mm protocol from WTCN
                            EchoSpacing = 540;
                        case 'b1epi2d3d2' % 800um protocol from WTCN
                            EchoSpacing = 540;
                        otherwise
                            fprintf(1,'WARNING: B1mapping version unknown, trying to base our guess on PixelBandwidth.\n');
                            PixelBandwidth = get_metadata_val(mstruc,'BandwidthPerPixelRO');
                            switch PixelBandwidth
                                case 2300
                                    EchoSpacing = 540e-3;
                                case 3600
                                    EchoSpacing = 330e-3;
                                case 3550 % Allegra data
                                    EchoSpacing = 330e-3;
                                otherwise
                                    fprintf(1,'Giving up: using default EchoSpacing value = 540 us.\n');
                                    EchoSpacing = 2*140+260; % us
                            end
                    end
                    measPElin = get_metadata_val(mstruc,'PELinesPAT');
                    cRes = 1;
                    parLocation{cRes} = 'HardCodedParameter';
                    parValue{cRes} = EchoSpacing * measPElin * 0.001; % ms
                else
                    if strcmp(valEPI,'EP')
                        fprintf(1,['WARNING: Sequence defined as EPI but BandwidthPerPixelPhaseEncode\n' ...
                            'not defined. No value returned.\n']);
                    else
                        fprintf(1,['WARNING: This might not be an EPI sequence. \n' ...
                            'Could not work ou EPI readout duration.\n']);
                    end
                end
            end
            
        case 'EchoSpacing' % [ms]
            % Siemens-specific
            % This information is easily retrievable from standard EPI
            % sequences, where the "BandwidthPerPixelPhaseEncoding" is defined.
            % In Customer-written sequences, this parameter might be missing in
            % the header...
            
            % first check whether BandwidthPerPixelPhaseEncode is defined
            [nFieldFound, fieldList] = find_field_name(mstruc, 'BandwidthPerPixelPhaseEncode','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
            % Keep only first value if many
            if nFieldFound
                cRes = 1;
                parLocation{cRes} = nam{1};
                epiROduration = 1/val{1}*1000;
                measPElin = get_metadata_val(mstruc,'PELinesPAT');
                try
                    parValue{cRes} = epiROduration/measPElin;
                catch
                    fprintf(1,'WARNING: Cannot retrieve the EchoSpacing for the current sequence.\n');
                end
            else
                fprintf(1,['WARNING: BandwidthPerPixelPhaseEncode not defined for the current sequence.\n' ...
                    'The echo spacing cannot be retrieved unambiguously.\n']);
            end
            
        case 'B1mapNominalFAValues' % [deg] for al_B1mapping - version dependent!!
            valSEQ = get_metadata_val(mstruc, 'SequenceName');
            valPROT = get_metadata_val(mstruc, 'ProtocolName');
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','insensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                switch lower(valSEQ)
                    case 'b1v2d3d2' % VD13 Prisma data
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                        parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                    case 'b1epi4a3d2' % VA35 Allegra data
                        % wip parameters are sorted as follows:
                        % alFree: [EddyCurrentDelay MixingTime NoRefAverages DurationPer5Deg BWT_SE/STE_factor NoDummyScans CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr CrusherAmplitude MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                        parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                    case 'b1epi2b3d2' % 1mm protocol from WTCN
                        % wip parameters are sorted as follows:
                        % alFree: [EddyCurrentDelay Tmixing (?) DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans]
                        parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                        parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                    case 'b1epi2d3d2' % 800um protocol from WTCN
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                        parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                    otherwise
                        fprintf(1,'B1mapping version unknown (%s/%s). Give up guessing FA values.\n', valSEQ, valPROT);
                end
                if ~isempty(parLocation)
                    nmeas = get_metadata_val(mstruc,'NumberOfMeasurements');
                    parValue{cRes} = parValue{cRes}(1:nmeas)*0.5;
                end
            end
            
        case 'RFSpoilingPhaseIncrement' % [Hz] defined in al_B1mapping and mtflash3d sequences - version dependent!!
            valSEQ = get_metadata_val(mstruc, 'SequenceName');
            valPROT = get_metadata_val(mstruc, 'ProtocolName');
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','insensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                index = 0;
                switch lower(valSEQ)
                    case 'b1v2d3d2' % VD13 Prisma
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 6;
                    case 'b1epi4a3d2' % VA35 Allegra data
                        % wip parameters are sorted as follows:
                        % alFree: [EddyCurrentDelay MixingTime NoRefAverages DurationPer5Deg BWT_SE/STE_factor NoDummyScans CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr CrusherAmplitude MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 6;
                    case 'fl3d_2l3d8' % VD13 Prisma
                        % wip parameters are sorted as follows:
                        % alFree: [RawDataExport(off/on=1/2) MTRepFactor DurationMTGaussianPulse FlatTopMTSpoiler DurPrewRamp DurPrewFlat DurRORamp RFExc(RectNonSel/SincNonSel/SincSlabSel = 1/2/3) RectFixedDur SincFixedDur BWTSinc]
                        % adFree: [MTGaussianFA OffResonanceMTGaussianPulse RFSpoilIncr]
                        index = 3;
                    case 'b1epi2d3d2' % 800um protocol from WTCN
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 6;                        
                    case 'fl3d_2d3d6' % VA35 Allegra
                        % wip parameters are sorted as follows:
                        % alFree: [MTSaturationMode (1/2 = Gaussian/Binomial)
                        %          MTRepFactor
                        %          BalancedMTSaturation (1/2 = false/true)
                        %          DurationMTGaussianPulse
                        %          RFExc(RectNonSel/SincNonSel/SincSlabSel = 1/2/3)
                        %          GRAPPA&RefScans (1/2 = false/true)
                        %          DurPrewRamp
                        %          DurPrewFlat
                        %          DurRORamp
                        %          FlatTopSpoiler]
                        % adFree: [MTGaussianFA OffResonanceMTGaussianPulse RFSpoilIncr SpoilerAmpl]
                        index = 3;
                    otherwise
                        fprintf(1,'Sequence version unknown (%s/%s). Give up guessing RF spoiling increment.\n', valSEQ, valPROT);
                end
                if index
                    parLocation{cRes} = [nam{1} '.adFree(' num2str(index) ')'];
                    try
                        parValue{cRes} = val{1}.adFree(index); % in deg
                    catch
                        % if index^th element = 0, it is not specified in the
                        % header and index exceeds matrix dimension:
                        parValue{cRes} = 0;
                    end
                end
            end
            
        case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
            valSEQ = get_metadata_val(mstruc, 'SequenceName');
            valPROT = get_metadata_val(mstruc, 'ProtocolName');
            [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','insensitive','matchType','exact');
            [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                index = 0;
                switch lower(valSEQ)
                    case 'b1v2d3d2'
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 1;
                    case 'b1epi4a3d2' % VA35 Allegra data
                        % wip parameters are sorted as follows:
                        % alFree: [EddyCurrentDelay MixingTime NoRefAverages DurationPer5Deg BWT_SE/STE_factor NoDummyScans CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr CrusherAmplitude MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 2;
                    case 'b1epi2b3d2' % 1mm protocol from WTCN
                        % wip parameters are sorted as follows:
                        % alFree: [EddyCurrentDelay Tmixing (?) DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans]
                        index = 2;
                    case 'b1epi2d3d2' % 800um protocol from WTCN
                        % wip parameters are sorted as follows:
                        % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                        % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                        index = 1;
                    otherwise
                        fprintf(1,'B1mapping version unknown (%s/%s). Give up guessing TM value.\n', valSEQ, valPROT);
                end
                if index
                    parLocation{cRes} = [nam{1} '.alFree(' num2str(index) ')'];
                    parValue{cRes} = val{1}.alFree(index)*0.001; % in ms
                end
            end
            
        otherwise
            [nFieldFound, fieldList] = find_field_name(mstruc, inParName, 'caseSens','insensitive','matchType','partial');
            [parValue,parLocation] = get_val_nam_list(mstruc, nFieldFound, fieldList);
            
    end
    if ~nFieldFound
        fprintf(1,'WARNING: No %s found in the extended header\n', inParName);
        parValue = [];
        parLocation = [];
    end
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
        cF = 1;
        nam{cRes} = fieldList{cRes,cF};
        while cF<size(fieldList,2)
            cF = cF+1;
            if ~isempty(fieldList{cRes,cF})
                nam{cRes} = [nam{cRes}  '.' fieldList{cRes,cF}];
            end
        end
        val{cRes} = eval(['mstruc.' nam{cRes}]);
    end
else
    val = {};
    nam = {};
end
end

