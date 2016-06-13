function [parValue, outParName] = hMRI_get_extended_hdr_val(hdr, inParName)
% PURPOSE
% To retrieve parameter values (mainly acquisition parameters, but this can
% be developed further in the future) from nifti extended header. The
% argument must be either a predefined parameter name (see examples below),
% or a potential match for a given field in the extended header.
%
% ARGUMENTS
% - hdr is the extended header structure
% - inParName can be an arbitrary string that will be searched in the
%   header structure or one of the following predefined parameter name:
%    - 'RepetitionTime' TR [ms]
%    - 'EchoTime' TE(s) [ms]
%    - 'FlipAngle' [deg]
%    - 'ProtocolName'
%    - 'SequenceName'
%    - 'MT' (1/0 = ON/OFF)
%    - 'FieldStrength' [T]
%    - 'Frequency' [Hz]
%    - 'ScanningSequence'
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
%   parName.
% - outParName is a string (or a cellarray of strings if the solution is
%   not unique) giving the location of the retrived value(s) in the header.
%
% SYNTAX
% [parValue, outParName] = hMRI_get_extended_hdr_val(hdr, 'RepetitionTime')
% [parValue, outParName] = hMRI_get_extended_hdr_val(hdr, 'Repetition')
% parValue = hMRI_get_extended_hdr_val(hdr, 'RepetitionTime')
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - June 2016 - Cyclotron Research Centre
%--------------------------------------------------------------------------

outParName = {};
nFieldFound = 0;

switch inParName
    case 'RepetitionTime' % [ms]
        [nFieldFound, fieldList] = findFieldName(hdr, 'alTR', 'caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Repetition time is given in us, we want it in ms
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1}*0.001;
        end
        
    case 'EchoTime' % [ms]
        [nFieldFound, fieldList] = findFieldName(hdr, 'alTE', 'caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Echo time is given in us, we want it in ms
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1}*0.001;
        end
        
    case 'FlipAngle' % [deg]
        [nFieldFound, fieldList] = findFieldName(hdr, 'FlipAngle', 'caseSens','sensitive','matchType','exact');
        %[nFieldFound, fieldList] = findFieldName(hdr, 'adFlipAngleDegree', 'caseSens','sensitive','matchType','exact'); % equivalent
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many - no scaling necessary
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'ProtocolName'
        [nFieldFound, fieldList] = findFieldName(hdr, 'ProtocolName','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'SequenceName'
        [nFieldFound, fieldList] = findFieldName(hdr, 'SequenceName','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'MT'
        % NB: parameters set to 0 are usually omitted in the DICOM header.
        % Therefore, the absence of the parameter in the header means that
        % no MT pulse is applied. If applied, parameter
        % acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPrepPulses.ucMTC
        % takes the hexadecimal value '0x1' = 1.
        [nFieldFound, fieldList] = findFieldName(hdr, 'ucMTC','caseSens','sensitive','matchType','exact');
        [~,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = 1;
        else
            nFieldFound = 1;
            cRes = 1;
            outParName{cRes} = 'NullParameterNotDefinedInStruct';
            parValue{cRes} = 0;
        end
        
    case 'FieldStrength' % [T]
        % NB: flNominalB0 returns ~2.8936 for a 3T magnet
        [nFieldFound, fieldList] = findFieldName(hdr, 'flNominalB0','caseSens','sensitive','matchType','exact');
        % while MagneticFieldStrength returns 3 for a 3T magnet
        % [nFieldFound, fieldList] = findFieldName(hdr, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'Frequency' % [Hz]
        % NB: lFrequency returns 123255074 Hz
        [nFieldFound, fieldList] = findFieldName(hdr, 'lFrequency','caseSens','sensitive','matchType','exact');
        % while ImagingFrequency returns 123.2551 MHz
        % [nFieldFound, fieldList] = findFieldName(hdr, 'MagneticFieldStrength','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'ScanningSequence' % e.g. 'EP' for EPI...
        [nFieldFound, fieldList] = findFieldName(hdr, 'ScanningSequence','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'MeasuredPELines' % taking PAT acceleration factor into account
        [nFieldFound, fieldList] = findFieldName(hdr, 'lPhaseEncodingLines','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        [nFieldFoundPAT, fieldListPAT] = findFieldName(hdr, 'lAccelFactPE','caseSens','sensitive','matchType','exact');
        [valPAT,~] = get_val_nam_list(hdr, nFieldFoundPAT, fieldListPAT);
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            if isempty(valPAT);valPAT{1} = 1;end
            parValue{cRes} = floor(val{1}/valPAT{1});
        end
        
    case 'PhaseEncodingDirectionSign'
        [nFieldFound, fieldList] = findFieldName(hdr, 'PhaseEncodingDirectionPositive','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        else
            nFieldFound = 1;
            cRes = 1;
            outParName{cRes} = 'NullParameterNotDefinedInStruct';
            parValue{cRes} = -1;
        end
        
    case 'PhaseEncodingDirection' % 'COL' (A>>P/P>>A) or 'ROW' (R>>L/L>>R)
        [nFieldFound, fieldList] = findFieldName(hdr, 'InPlanePhaseEncodingDirection','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'NumberOfMeasurements'
        [nFieldFound, fieldList] = findFieldName(hdr, 'lRepetitions','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1}+1;
        else
            nFieldFound = 1;
            cRes = 1;
            outParName{cRes} = 'NullParameterNotDefinedInStruct';
            parValue{cRes} = 1;
        end
        
    case 'WipParameters'
        [nFieldFound, fieldList] = findFieldName(hdr, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = val{1};
        end
        
    case 'AllDiffusionDirections'
        if hMRI_get_extended_hdr_val(hdr,'isDWI')
            [nFieldFound, fieldList] = findFieldName(hdr, 'sDiffusion','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
            % sDiffusion is the field containing series diffusion information.
            % Example:
            % hdr{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sDiffusion
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
            % hdr{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sDiffusion.sFreeDiffusionData
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
                    ndir = eval(['hdr.' nam{1} '.sFreeDiffusionData.lDiffDirections']);
                    outParName{cRes} = [nam{1} '.sFreeDiffusionData.asDiffDirVector'];
                    parValueCell = eval(['hdr.' nam{1} '.sFreeDiffusionData.asDiffDirVector']);
                    parValue{cRes} = zeros(3,ndir);
                    for cdir = 1:length(parValueCell)
                        if isempty(parValueCell{cdir}.dSag);parValueCell{cdir}.dSag = 0;end
                        if isempty(parValueCell{cdir}.dCor);parValueCell{cdir}.dCor = 0;end
                        if isempty(parValueCell{cdir}.dTra);parValueCell{cdir}.dTra = 0;end
                        parValue{cRes}(:,cdir) = [parValueCell{cdir}.dSag; parValueCell{cdir}.dCor; parValueCell{cdir}.dTra];
                    end
                end
                
            end
        end
        
    case 'AllBValues'
        if hMRI_get_extended_hdr_val(hdr,'isDWI')
            [nFieldFound, fieldList] = findFieldName(hdr, 'alBValue','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
            
            
            if nFieldFound
                cRes = 1;
                outParName{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
        end
        
    case 'DiffusionDirection'
        if hMRI_get_extended_hdr_val(hdr,'isDWI')
            [nFieldFound, fieldList] = findFieldName(hdr, 'DiffusionGradientDirection','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                outParName{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
        end
        
    case 'BValue'
        if hMRI_get_extended_hdr_val(hdr,'isDWI')
            [nFieldFound, fieldList] = findFieldName(hdr, 'B_value','caseSens','sensitive','matchType','exact');
            [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
            if nFieldFound
                cRes = 1;
                outParName{cRes} = nam{1};
                parValue{cRes} = val{1};
            end
        end
        
    case 'isDWI'
        [nFieldFound, fieldList] = findFieldName(hdr, 'sDiffusion','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            outParName{cRes} = [nam{1} '.ulMode'];
            if (val{1}.ulMode>1)
                parValue{cRes} = 1;
            else
                parValue{cRes} = 0;
            end
        else
            nFieldFound = 1;
            cRes = 1;
            outParName{cRes} = 'NotDWISequence';
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
        [nFieldFound, fieldList] = findFieldName(hdr, 'BandwidthPerPixelPhaseEncode','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        % if (nFieldFound>1);warning('More than one value was found for %s. First one kept.', inParName);end
        % Keep only first value if many
        if nFieldFound
            cRes = 1;
            outParName{cRes} = nam{1};
            parValue{cRes} = 1/val{1}*1000;
        else
            warning(['BandwidthPerPixelPhaseEncode not defined for the current sequence\n' ...
                'For 3D-EPI B1 mapping sequences, values are deduced from the sequence\n' ...
                'version. Be aware that it might not be correct if the version is unknown']);
            % check whether it is an EPI sequence (al_B1mapping is not defined as 'EP' but 'RM'):
            valEPI = hMRI_get_extended_hdr_val(hdr, 'ScanningSequence');
            valSEQ = hMRI_get_extended_hdr_val(hdr, 'SequenceName');
            valPROT = hMRI_get_extended_hdr_val(hdr, 'ProtocolName');
            if strcmp(valEPI,'EP')
                warning('Sequence defined as EPI but BandwidthPerPixelPhaseEncode. No value returned.');
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
                measPElin = hMRI_get_extended_hdr_val(hdr,'MeasuredPELines');
                cRes = 1;
                outParName{cRes} = 'HardCodedParameter';
                parValue{cRes} = EchoSpacing * measPElin * 0.001; % ms
            end
        end
        
    case 'B1mapNominalFAValues' % [deg] for al_B1mapping - version dependent!!
        valSEQ = hMRI_get_extended_hdr_val(hdr, 'SequenceName');
        valPROT = hMRI_get_extended_hdr_val(hdr, 'ProtocolName');
        [nFieldFound, fieldList] = findFieldName(hdr, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    outParName{cRes} = [nam{1} '.adFree(3:4)'];
                    parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                otherwise
                    warning('B1mapping version unknown (%s / %s). Give up guessing FA values.', valSEQ, valPROT);
            end
            if ~isempty(outParName)
                nmeas = hMRI_get_extended_hdr_val(hdr,'NumberOfMeasurements');
                parValue{cRes} = parValue{cRes}(1:nmeas)*0.5;
            end
        end
        
    case 'RFSpoilingPhaseIncrement' % [Hz] defined in al_B1mapping and mtflash3d sequences - version dependent!!
        valSEQ = hMRI_get_extended_hdr_val(hdr, 'SequenceName');
        valPROT = hMRI_get_extended_hdr_val(hdr, 'ProtocolName');
        [nFieldFound, fieldList] = findFieldName(hdr, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    outParName{cRes} = [nam{1} '.adFree(6)'];
                    parValue{cRes} = val{1}.adFree(6); % in deg
                case 'fl3d_2l3d8'
                    % wip parameters are sorted as follows:
                    % alFree: [RawDataExport(off/on=1/2) MTRepFactor DurationMTGaussianPulse FlatTopMTSpoiler DurPrewRamp DurPrewFlat DurRORamp RFExc(RectNonSel/SincNonSel/SincSlabSel = 1/2/3) RectFixedDur SincFixedDur BWTSinc]
                    % adFree: [MTGaussianFA OffResonanceMTGaussianPulse RFSpoilIncr]
                    outParName{cRes} = [nam{1} '.adFree(3)'];
                    parValue{cRes} = val{1}.adFree(3); % in deg
                otherwise
                    warning('Sequence version unknown (%s / %s). Give up guessing RF spoiling increment.', valSEQ, valPROT);
            end
        end
        
    case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
        valSEQ = hMRI_get_extended_hdr_val(hdr, 'SequenceName');
        valPROT = hMRI_get_extended_hdr_val(hdr, 'ProtocolName');
        [nFieldFound, fieldList] = findFieldName(hdr, 'sWipMemBlock','caseSens','sensitive','matchType','exact');
        [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList);
        if nFieldFound
            cRes = 1;
            switch valSEQ
                case 'B1v2d3d2'
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor CursherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    outParName{cRes} = [nam{1} '.alFree(1)'];
                    parValue{cRes} = val{1}.alFree(1)*0.001; % in ms
                otherwise
                    warning('B1mapping version unknown (%s / %s). Give up guessing TM value.', valSEQ, valPROT);
            end
        end
        
    otherwise
        [nFieldFound, fieldList] = findFieldName(hdr, inParName, 'caseSens','insensitive','matchType','partial');
        [parValue,outParName] = get_val_nam_list(hdr, nFieldFound, fieldList);
        
end

if ~nFieldFound
    warning('No %s found in the extended header', inParName);
    parValue = [];
    outParName = {};
end

% returns cell array only if necessary (non-unique result)
if length(parValue) == 1
    parValue = parValue{1};
    outParName = outParName{1};
end
end


function [val,nam] = get_val_nam_list(hdr, nFieldFound, fieldList) %#ok<INUSL>
% to retrieve the values and location in the hdr structure
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
        val{cRes} = eval(['hdr.' nam{cRes}]);
    end
else
    val = {};
    nam = {};
end
end

