function [value, result] = get_val_classic(obj, fieldName)
  value = [];
  result = false;

  switch fieldName
    case 'MT'
      % Siemens-specific
      % NB: parameters set to 0 are usually omitted in the DICOM header.
      % Therefore, the absence of the parameter in the header means that
      % no MT pulse is applied. If applied, parameter
      % acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPrepPulses.ucMTC
      % takes the hexadecimal value '0x1' = 1.
      [value, result] = obj.get_val_raw(obj.meta_hash, 'ucMTC', true);
      if result
        value = double(value > 0);
      else
        reuslt = true;
        value = 0;
      end
        
    case 'FieldStrength' % [T]
      [value, result] = obj.get_val_raw(obj.meta_hash, 'flNominalB0', true);
      if ~result
        [value, result] = obj.get_val_raw(obj.meta_hash, 'MagneticFieldStrength', true);
      end
        
    case 'Frequency' % [Hz]
      % Valid for all vendors
      % NB: lFrequency returns 123255074 Hz and is Siemens-specific,
      % while ImagingFrequency returns 123.2551 MHz.
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lFrequency', true);
      if ~result
        [value, result] = obj.get_val_raw(obj.meta_hash, 'MagneticFieldStrength', true);
        if result
          value = value * 1e6;
        end
      end
        
    case 'BandwidthPerPixelRO'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'PixelBandwidth', true);
        
    case 'BandwidthPerPixelPE' % Siemens-specific header entry
      [value, result] = obj.get_val_raw(obj.meta_hash, 'BandwidthPerPixelPhaseEncode', true);
        
    case 'PELinesPF'
      % size of the k-space PE dimension, taking into account partial
      % Fourier but not Parallel acceleration. Valid for all vendors.
      [value, result] = obj.get_val_raw(obj.meta_hash, 'NumberOfPhaseEncodingSteps', true);
        
    case 'PELines'
      % size of the k-space PE dimension, without taking into account
      % Parallel acceleration nor partial Fourier.
      % Siemens-specific towards validation for all vendors...
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lPhaseEncodingLines', true);
      if ~result
        [value, result] = obj.get_val_classic('PhaseEncodingDirection');
        if strcmp(value, 'ROW')
          [value, result] = obj.get_val_raw(obj.meta_hash, 'Rows');
        else
          [value, result] = obj.get_val_raw(obj.meta_hash, 'Columns');
        end
      end
        
    case 'PELinesPAT'
      % size of the k-space PE dimension, taking into account Parallel
      % acceleration but not partial Fourier. Used to calculate the total
      % EPI Readout duration for FieldMap undistortion.
      % Siemens-specific.
      [value, result] = obj.get_val_classic('PELines', true);
      [valuePAT, resultPAT] = obj.get_val_raw(obj.meta_hash, 'lAccelFactPE', true);
      if ~resultPAT
        result = true;
        value = 1;
      elseif result
        value = floor(value / valuePAT);
      end
        
    case 'PhaseEncodingDirectionSign'
      % Siemens-specific:
      [value, result] = obj.get_val_raw(obj.meta_hash, 'PhaseEncodingDirectionPositive', true);
      if ~result
        result = true;
        value = -1;
      end
        
    case 'NumberOfMeasurements' % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lRepetitions', true);
      if ~result
        result = true;
        value = 1;
      end
        
    case 'NumberOfSlices' % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'sSliceArray', true);
      if result
        if value.lSize == 1 && ...
            strcmp(obj.get_val_raw(obj.meta_hash, 'MRAcquisitionType'),'3D')
          [value, result] = obj.get_val_raw(obj.meta_hash, 'lImagesPerSlab', true);
        else
          value = value.lSize;
        end
      end
        
    case 'PATparameters'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'sPat', true);
        
    case 'AccelFactorPE'
      % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lAccelFactPE', true);
        
    case 'AccelFactor3D'
      % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lAccelFact3D', true);
        
    case 'MultiBandFactor'
      % Siemens-specific for either CMRR multiband EPI or Siemens
      % product EPI sequence with SMS.
      % determine sequence (CMRR versus Siemens product)
      [valSeq, resSeq] = obj.get_val_raw(obj.meta_hash, 'tSequenceFileName', true);
      if resSeq
        if strfind(lower(valSeq), 'cmrr')
          [value, result] = obj.get_val_raw(obj.meta_hash, 'MiscSequenceParam', true);
          if result
            value = value(12);
          end
        else
          [value, result] = obj.get_val_raw(obj.meta_hash, 'lMultiBandFactor', true);
        end
      end 
        
    case 'WipParameters'
      % Siemens-specific (NB: search made case insensitive since
      % sWiPMemBlock or sWipMemBlock depending on software version)
      [value, result] = obj.get_val_raw(obj.meta_hash, 'sWipMemBlock', false);
        
    case 'DiffusionDirection'
      % Siemens-specific
      if obj.get_val_classic('isDWI')
        [value, result] = obj.get_val_raw(obj.meta_hash, 'DiffusionGradientDirection', false);
        if ~result
          value = [0; 0; 0];
          result = true;
        end
      end
        
    case 'BValue'
        % Siemens-specific
        if obj.get_val_classic('isDWI')
          [value, result] = obj.get_val_raw(obj.meta_hash, 'B_value', false);
        end
        
    case 'isDWI'
        % Siemens-specific
        [value, result] = obj.get_val_raw(obj.meta_hash, 'sDiffusion', false);
        if result
          value = double(value.ulMode > 1);
        else
          result = true;
          value = 0;
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
            fprintf(1,['\nWARNING: BandwidthPerPixelPhaseEncode not defined for the current sequence\n' ...
                'For 3D-EPI B1 mapping sequences, values are deduced from the sequence\n' ...
                'version. Be aware that it might not be correct if the version is unknown.\n']);
            
            % check whether it is an EPI sequence (al_B1mapping is not defined as 'EP' but 'RM'):
            valEPI = get_metadata_val(mstruc, 'ScanningSequence');
            valSEQ = get_metadata_val(mstruc, 'SequenceName');
            valPROT = get_metadata_val(mstruc, 'ProtocolName');
            valMODELNAME = get_metadata_val(mstruc, 'ManufacturerModelName');
            
            % Case al_B1mapping
            if strfind(lower(valPROT),'b1map')
                fprintf(1,'\nTrying to derive the epiReadoutDuration from the Sequence version (%s/%s).\n',valSEQ,valPROT);
                nFieldFound = 1;
                switch lower(valSEQ)
                    case 'b1v2d3d2' % 540 us - Prisma
                        EchoSpacing = 2*140+260;
                    case 'b1epi4a3d2' % 330 us - Allegra
                        EchoSpacing = 330;
                    case 'b1epi2b3d2' % 1mm protocol from WTCN
                        EchoSpacing = 540;
                    case 'seste1d3d2' % 1mm protocol from WTCN
                        % alFree: [VoxDeph,SpoilAmp,EddCurr0,EddCurr1,TRamp,TFlat,BWT,0,0,0,0,0,2,MixingTime,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12345],
                        % adFree: [0,0,0,0,0,0,0,SlabGradScale,RefocCorr,0,0,RFSpoilBasicIncr]
                        Wip = get_metadata_val(mstruc, 'WipParameters');
                        EchoSpacing = Wip.alFree(5)*2+Wip.alFree(6);
                    case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'} % 7T and Prisma versions by Kerrin Pine
                        if contains(valMODELNAME,'Prisma','IgnoreCase',true)
                            Wip = get_metadata_val(mstruc, 'WipParameters');
                            EchoSpacing = Wip.alFree(5)*2+Wip.alFree(6);
                        elseif contains(valMODELNAME,'7T','IgnoreCase',true)
                            EchoSpacing = 540;
                        else
                            % Do nothing
                        end
                    case 'b1epi2d3d2' % 800um protocol from WTCN
                        EchoSpacing = 540;
                end
                if ~exist('EchoSpacing','var')
                    fprintf(1,'\nWARNING: B1mapping version unknown, trying to base our guess on PixelBandwidth.\n');
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
                
                fprintf(1,['\nINFO: the EPI readout duration has been derived:' ...
                    '\n\tEchoSpacing = %5.2f us' ...
                    '\n\tmeasPElin = %d' ...
                    '\n\tepiReadoutDuration = %5.2f ms\n'], ...
                    EchoSpacing, measPElin,parValue{cRes});
                
            else
                if strcmp(valEPI,'EP')
                    fprintf(1,['\nWARNING: Sequence defined as EPI but BandwidthPerPixelPhaseEncode\n' ...
                        'not defined. No value returned.\n']);
                else
                    fprintf(1,['\nWARNING: This might not be an EPI sequence. \n' ...
                        'Could not work out EPI readout duration.\n']);
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
                fprintf(1,'\nWARNING: Cannot retrieve the EchoSpacing for the current sequence.\n');
            end
        else
            fprintf(1,['\nWARNING: BandwidthPerPixelPhaseEncode not defined for the current sequence.\n' ...
                'Trying to derive the echo spacing from the epiReadoutDuration...\n']);
            cRes = 1;
            epiROduration = get_metadata_val(mstruc,'epiReadoutDuration');
            measPElin = get_metadata_val(mstruc,'PELinesPAT');
            try
                parValue{cRes} = epiROduration/measPElin;
                parLocation{cRes} = 'Derived';
                nFieldFound = 1;
            catch
                fprintf(1,'\nWARNING: Failed deriving the Echo Spacing for the current protocol.\n');
            end
        end
        
    case 'B1mapNominalFAValues' % [deg] for al_B1mapping - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        valMODELNAME = get_metadata_val(mstruc, 'ManufacturerModelName');
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
                case 'seste1d3d2' % 1mm protocol from WTCN (MFC)
                    % wip parameters are sorted as follows:
                    % alFree: [VoxDeph,SpoilAmp,EddCurr0,EddCurr1,TRamp,TFlat,BWT,0,0,0,0,0,2,MixingTime,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12345],
                    % adFree: [0,0,0,0,0,0,0,SlabGradScale,RefocCorr,0,0,RFSpoilBasicIncr]
                    parLocation{cRes} = 'HardCodedParameter';
                    parValue{cRes} = 230:-10:0;
                case 'b1epi2d3d2' % 800um protocol from WTCN
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                    parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;
                case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'} % 7T and Prisma versions by Kerrin Pine
                    if contains(valMODELNAME,'Prisma','IgnoreCase',true)
                        parLocation{cRes} = 'HardCodedParameter';
                        parValue{cRes} = 230:-10:0;
                    elseif contains(valMODELNAME,'7T','IgnoreCase',true)
                        parLocation{cRes} = [nam{1} '.adFree(3:4)'];
                        parValue{cRes} = val{1}.adFree(3):-val{1}.adFree(4):0;      
                    else
                        parLocation{cRes}=[];
                    end
                otherwise
                    fprintf(1,'\nWARNING: B1mapping version unknown (%s/%s). Give up guessing FA values.\n', valSEQ, valPROT);
            end
            if ~isempty(parLocation)
                nmeas = get_metadata_val(mstruc,'NumberOfMeasurements');
                parValue{cRes} = parValue{cRes}(1:nmeas)*0.5;
            end
        end
        
    case 'RFSpoilingPhaseIncrement' % [deg] defined in al_B1mapping and mtflash3d sequences - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        valMODELNAME = get_metadata_val(mstruc, 'ManufacturerModelName');
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
                case 'seste1d3d2' % 1mm protocol from WTCN (MFC)
                    % wip parameters are sorted as follows:
                    % alFree: [VoxDeph,SpoilAmp,EddCurr0,EddCurr1,TRamp,TFlat,BWT,0,0,0,0,0,2,MixingTime,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12345],
                    % adFree: [0,0,0,0,0,0,0,SlabGradScale,RefocCorr,0,0,RFSpoilBasicIncr]
                    index = 12;
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
                case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'} % Prisma and 7T versions by Kerrin Pine
                    if contains(valMODELNAME,'Prisma','IgnoreCase',true)
                        index = 12;
                    elseif contains(valMODELNAME,'7T','IgnoreCase',true)
                        index = 3;
                    else
                        index = false;
                    end
                otherwise
                    fprintf(1,'Sequence version unknown (%s/%s). Give up guessing RF spoiling increment.\n', valSEQ, valPROT);
            end
            if index
                parLocation{cRes} = [nam{1} '.adFree(' num2str(index) ')'];
                try
                    parValue{cRes} = val{1}.adFree(index); % in deg
                catch %#ok<*CTCH>
                    % if index^th element = 0, it is not specified in the
                    % header and index exceeds matrix dimension:
                    parValue{cRes} = 0;
                end
            end
        end
        
    case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
        valSEQ = get_metadata_val(mstruc, 'SequenceName');
        valPROT = get_metadata_val(mstruc, 'ProtocolName');
        valMODELNAME = get_metadata_val(mstruc, 'ManufacturerModelName');
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
                case 'seste1d3d2' % 1mm protocol from WTCN (MFC)
                    % wip parameters are sorted as follows:
                    % alFree: [VoxDeph,SpoilAmp,EddCurr0,EddCurr1,TRamp,TFlat,BWT,0,0,0,0,0,2,MixingTime,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12345],
                    % adFree: [0,0,0,0,0,0,0,SlabGradScale,RefocCorr,0,0,RFSpoilBasicIncr]
                    index = 14;
                case 'b1epi2d3d2' % 800um protocol from WTCN
                    % wip parameters are sorted as follows:
                    % alFree: [Tmixing DurationPer5Deg BWT_SE/STE_factor (?) CrusherPerm(on/off=2/3) OptimizedRFDur(on/off=2/3)]
                    % adFree: [RefocCorr ScaleSGrad? MaxRefocAngle DecRefocAngle FAforReferScans RFSpoilIncr]
                    index = 1;
                case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'} % 7T and Prisma versions by Kerrin Pine
                    if contains(valMODELNAME,'Prisma','IgnoreCase',true)
                        index = 14;
                    elseif contains(valMODELNAME,'7T','IgnoreCase',true)
                        index = 1;
                    else
                        index = false;
                    end
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
end
