function [parValue, parLocation] = get_metadata_val(varargin)
% This is hmri_get_metadata_val, part of the hMRI-Toolbox
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
%   metadata structure or one of the following predefined parameter name:
%    - 'RepetitionTime' TR [ms]
%    - 'EchoTime' TE(s) [ms]
%    - 'FlipAngle' [deg]
%    - 'ProtocolName'
%    - 'SequenceName'
%    - 'MT' (1/0 = ON/OFF)
%    - 'FieldStrength' [T]
%    - 'Frequency' [Hz]
%    - 'ScanningSequence'
%    - 'BandwidthPerPixelRO' [Hz/Px]
%    - 'BandwidthPerPixelPE' [Hz/Px]
%    - 'PATparameters' [struct]
%    - 'AccelFactorPE'
%    - 'AccelFactor3D'
%    - 'MeasuredPELines'
%    - 'PhaseEncodingDirectionPositive' A>>P & R>>L = 1; P>>A & L>>R = 0.
%    - 'PhaseEncodingDirection' A>>P/P>>A = 'COL' & R>>L/L>>R = 'ROW'
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
%   cellarray, ...). If parName is not leading to a unique value, a cell
%   array of values is returned, each element corresponding to a different
%   location in the header structure, as specified by the returned
%   parLocation.
% - parLocation is a string (or a cellarray of strings if the solution is
%   not unique) giving the location of the retrieved value(s) in the
%   metadata structure.
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

switch inParName
    case 'RepetitionTime' % [ms]
        [nFieldFound, fieldList] = find_field_name(mstruc, 'alTR', 'caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Repetition time is given in us, we want it in ms
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1}*0.001;
        end
        
    case 'EchoTime' % [ms]
        [nFieldFound, fieldList] = find_field_name(mstruc, 'alTE', 'caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Echo time is given in us, we want it in ms
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1}*0.001;
        end
        
    case 'FlipAngle' % [deg]
        [nFieldFound, fieldList] = find_field_name(mstruc, 'FlipAngle', 'caseSens','sensitive','matchType','exact');
        %[nFieldFound, fieldList] = find_field_name(mstruc, 'adFlipAngleDegree', 'caseSens','sensitive','matchType','exact'); % equivalent
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many - no scaling necessary
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'ProtocolName'
        [nFieldFound, fieldList] = find_field_name(mstruc, 'ProtocolName','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'SequenceName'
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
        % NB: flNominalB0 returns ~2.8936 for a 3T magnet
        [nFieldFound, fieldList] = find_field_name(mstruc, 'flNominalB0','caseSens','sensitive','matchType','exact');
        % while MagneticFieldStrength returns 3 for a 3T magnet
        % [nFieldFound, fieldList] = find_field_name(mstruc, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'Frequency' % [Hz]
        % NB: lFrequency returns 123255074 Hz
        [nFieldFound, fieldList] = find_field_name(mstruc, 'lFrequency','caseSens','sensitive','matchType','exact');
        % while ImagingFrequency returns 123.2551 MHz
        % [nFieldFound, fieldList] = find_field_name(mstruc, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'ScanningSequence' % e.g. 'EP' for EPI...
        [nFieldFound, fieldList] = find_field_name(mstruc, 'ScanningSequence','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'BandwidthPerPixelRO' % e.g. 'EP' for EPI...
        [nFieldFound, fieldList] = find_field_name(mstruc, 'PixelBandwidth','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'BandwidthPerPixelPE' % e.g. 'EP' for EPI...
        [nFieldFound, fieldList] = find_field_name(mstruc, 'BandwidthPerPixelPhaseEncode','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'MeasuredPELines' % taking PAT acceleration factor into account
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
        [nFieldFound, fieldList] = find_field_name(mstruc, 'PhaseEncodingDirectionPositive','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
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
        [nFieldFound, fieldList] = find_field_name(mstruc, 'InPlanePhaseEncodingDirection','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'NumberOfMeasurements'
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
        [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            parLocation{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'AllDiffusionDirections'
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
                warning('Diffusion direction not defined for DWImage %s. Assuming b=0.', inParName);
            end
        end
        
    case 'BValue'
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
            warning('This is not a DWI sequence');
        end
        
        
    case 'epiReadoutDuration' % [ms]
        % This information is easily retrievable from standard EPI
        % sequences, where the "BandwidthPerPixelPhaseEncoding" is defined.
        % For the 3D-EPI B1mapping sequence from Antoine, everything is
        % hard coded in the sequence and not passed to the MrProt variable,
        % therefore it is not available in the DICOM header. We have to
        % work case by case, according to the sequence version and
        % resolution :/...
        
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
            warning(['BandwidthPerPixelPhaseEncode not defined for the current sequence\n' ...
                'For 3D-EPI B1 mapping sequences, values are deduced from the sequence\n' ...
                'version. Be aware that it might not be correct if the version is unknown']);
            % check whether it is an EPI sequence (al_B1mapping is not defined as 'EP' but 'RM'):
            valEPI = get_metadata_val(mstruc, 'ScanningSequence');
            valSEQ = get_metadata_val(mstruc, 'SequenceName');
            valPROT = get_metadata_val(mstruc, 'ProtocolName');
            if strcmp(valEPI,'EP')
                warning('Sequence defined as EPI but BandwidthPerPixelPhaseEncode not defined. No value returned.');
            elseif strfind(lower(valSEQ),'b1')
                warning('Trying to derive the epiReadoutDuration from the Sequence version (%s / %s)',valSEQ,valPROT);
                nFieldFound = 1;
                switch valSEQ
                    case 'B1v2d3d2'
                        EchoSpacing = 2*140+260; % us
                    otherwise
                        warning('B1mapping version unknown, using default EchoSpacing value = 540 us.')
                        EchoSpacing = 2*140+260; % us
                end
                measPElin = get_metadata_val(mstruc,'MeasuredPELines');
                cRes = 1;
                parLocation{cRes} = 'HardCodedParameter';
                parValue{cRes} = EchoSpacing * measPElin * 0.001; % ms
            end
        end
       
        
    case 'EchoSpacing' % [ms]
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
            measPElin = get_metadata_val(mstruc,'MeasuredPELines');
            try
                parValue{cRes} = epiROduration/measPElin;
            catch
                warning('Cannot retrieve the EchoSpacing for the current sequence');
            end    
        else
            warning(['BandwidthPerPixelPhaseEncode not defined for the current sequence\n' ...
                'The echo spacing cannot be retrieved unambiguously.']);
        end
        
    case 'B1mapNominalFAValues' % [deg] for al_B1mapping - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                    parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                otherwise
                    warning('B1mapping version unknown (%s / %s). Give up guessing FA values.', valSEQ, valPROT);
            end
            if ~isempty(parLocation)
                nmeas = get_metadata_val(mstruc,'NumberOfMeasurements');
                parValue{cRes} = parValue{cRes}(1:nmeas)*0.5;
            end
        end
        
    case 'RFSpoilingPhaseIncrement' % [Hz] defined in al_B1mapping and mtflash3d sequences - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    parLocation{cRes} = [nam{1} '.adFree(6)'];
                    parValue{cRes} = val{1}.adFree(6); % in deg
                case 'fl3d_2l3d8'
                    % wip parameters are sorted as follows:
                    % alFree: [RawDataExport(off/on=1/2) MTRepFactor DurationMTGaussianPulse FlatTopMTSpoiler DurPrewRamp DurPrewFlat DurRORamp RFExc(RectNonSel/SincNonSel/SincSlabSel = 1/2/3) RectFixedDur SincFixedDur BWTSinc]
                    % adFree: [MTGaussianFA OffResonanceMTGaussianPulse RFSpoilIncr]
                    parLocation{cRes} = [nam{1} '.adFree(3)'];
                    parValue{cRes} = val{1}.adFree(3); % in deg
                otherwise
                    warning('Sequence version unknown (%s / %s). Give up guessing RF spoiling increment.', valSEQ, valPROT);
            end
        end
        
    case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        [nFieldFound, fieldList] = find_field_name(mstruc, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    parLocation{cRes} = [nam{1} '.alFree(1)'];
                    parValue{cRes} = val{1}.alFree(1)*0.001; % in ms
                otherwise
                    warning('B1mapping version unknown (%s / %s). Give up guessing TM value.', valSEQ, valPROT);
            end
        end
        
    otherwise
        [nFieldFound, fieldList] = find_field_name(mstruc, inParName, 'caseSens','insensitive','matchType','partial');
        [parValue,parLocation] = get_val_nam_list(mstruc, nFieldFound, fieldList);
        
end

if ~nFieldFound
    warning('No %s found in the extended header', inParName);
    parValue = [];
    parLocation = [];
end

% returns cell array only if necessary (non-unique result)
if length(parValue) == 1
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

