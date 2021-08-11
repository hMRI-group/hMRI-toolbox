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
        result = true;
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
        [valPE, resPE] = obj.get_val_classic('PhaseEncodingDirection');
        if resPE && strcmp(valPE, 'ROW')
          [value, result] = obj.get_val_raw(obj.meta_hash, 'Rows', true);
        else
          [value, result] = obj.get_val_raw(obj.meta_hash, 'Columns', true);
        end
      end
        
    case 'PELinesPAT'
      % size of the k-space PE dimension, taking into account Parallel
      % acceleration but not partial Fourier. Used to calculate the total
      % EPI Readout duration for FieldMap undistortion.
      % Siemens-specific.
      [value, result] = obj.get_val_classic('PELines');
      [valuePAT, resultPAT] = obj.get_val_raw(obj.meta_hash, 'lAccelFactPE', true);
      if ~resultPAT
        result = true;
        value = 1;
      elseif result
        value = floor(value / valuePAT);
      end
        
    case 'PhaseEncodingDirection'
      % Siemens-specific:
      [value, result] = obj.get_val_raw(obj.meta_hash, 'InPlanePhaseEncodingDirection', true);

    case 'PhaseEncodingDirectionSign'
      % Siemens-specific:
      [value, result] = obj.get_val_raw(obj.meta_hash, 'PhaseEncodingDirectionPositive', true);
      if ~result
        result = true;
        value = -1;
      end
        
    case 'NumberOfMeasurements' % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'lRepetitions', true);
      if result
        value = value + 1;
      else
        result = true;
        value = 1;
      end
        
    case 'NumberOfSlices' % Siemens-specific
      [value, result] = obj.get_val_raw(obj.meta_hash, 'sSliceArray', true);
      if result
        if value.lSize == 1 && ...
            strcmp(obj.get_val_raw(obj.meta_hash, 'MRAcquisitionType', true),'3D')
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
        if contains(lower(valSeq), 'cmrr')
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
        [value, result] = obj.get_val_raw(obj.meta_hash, 'BandwidthPerPixelPhaseEncode', true);
        if result
          % BPPPE is defined, deriving readout duration directly
          value = 1 / value * 1e3;
        else
          % BPPPE not defined, deriving readout duration from protocol
          valEPI = obj.get_val_raw(obj.meta_hash, 'ScanningSequence', true);
          valPROT = obj.get_val_raw(obj.meta_hash, 'ProtocolName', true);

          if contains(lower(valPROT),'b1map')
            [EchoSpacing, statusES] = obj.get_val_classic('EchoSpacing');
            [measPElin, statusPE] = obj.get_val_classic('PELinesPAT');
            if statusES && statusPE
              value = EchoSpacing * measPElin;
              result = true;
            end
          else
            if strcmp(valEPI,'EP')
              warning(['%s (%s): Sequence is EPI but ' ...
                       'BandwidthPerPixelPhaseEncode not defined'], ...
                      obj.filelist{obj.index}, fieldName);
            else
              warning('%s (%s): Sequence not EPI', ...
                      obj.filelist{obj.index}, fieldName);
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
      [BPPPE, resBPP] = obj.get_val_raw(obj.meta_hash, 'BandwidthPerPixelPhaseEncode', true);
      [PELines, resPELines] = obj.get_val_classic('PELinesPAT');
      if resBPP && resPELines
        value = 1. / (BPPPE * PELines);
        result = true;
      end

      if ~result
        % Retriving EchoSpacing based on Sequence
        valSEQ = obj.get_val_raw(obj.meta_hash, 'SequenceName', true);
        valMODELNAME = obj.get_val_raw(obj.meta_hash, 'ManufacturerModelName', true);
        wip = obj.get_val_raw(obj.meta_hash, 'sWipMemBlock', false);
        [value, result] = get_EchoSpacing(valSEQ, wip, valMODELNAME);
      end

      if ~result
        % Based on BPPPRO
        [PixelBandwidth, resPB]  = obj.get_val_raw(obj.meta_hash, 'BandwidthPerPixelRO', true);
        if resPB
          [value, result] = get_EchoSpacing_BPP(PixelBandwidth);
        end
      end

      if ~result
        % Default value
        value = 0.540; % 2*140+260 us
        result = true;
        warning('%s (%s): Unable to get EchoSpacing. Using default %f', ...
                obj.filelist{obj.index}, fieldName, value);
      end
        
    case 'RFSpoilingPhaseIncrement'
      % [deg] defined in al_B1mapping and mtflash3d sequences - version dependent!!
      valSEQ = obj.get_val_raw(obj.meta_hash, 'SequenceName', true);
      valPROT = obj.get_val_raw(obj.meta_hash, 'ProtocolName', true);
      valMODELNAME = obj.get_val_raw(obj.meta_hash, 'ManufacturerModelName', true);
      [wip, statusWIP] = obj.get_val_raw(obj.meta_hash, 'sWipMemBlock', false);
      index = 0;

      switch lower(valSEQ)
        case {'b1v2d3d2', 'b1epi4a3d2', 'b1epi2d3d2'} 
          index = 6;

        case {'fl3d_2l3d8', 'fl3d_2d3d6'}
          index = 3;

        case {'seste1d3d2'}
          index = 12;

        case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'}
          if contains(valMODELNAME,'Prisma','IgnoreCase',true)
              index = 12;
          elseif contains(valMODELNAME,'7T','IgnoreCase',true)
              index = 3;
          end
        otherwise
          warning(['%s (%s): Unknown sequence %s (%s). ' ...
                   'Unable to extract RF spoiling increment'], ...
                  obj.filelist{obj.index}, fieldName, valSEQ, valPROT)
      end
      if statusWIP && index > 0
        if index <= numel(wip.adFree)
          value = wip.adFree(index);
          result = true;
        else
          warning('%s (%s): Index %d out of range for sWipMemBlock', ...
                  obj.filelist{obj.index}, fieldName, index);
        end
      end
        
    case 'B1mapMixingTime' % [ms] for al_B1mapping - version dependent!!
      valSEQ = obj.get_val_raw(obj.meta_hash, 'SequenceName', true);
      valPROT = obj.get_val_raw(obj.meta_hash, 'ProtocolName', true);
      valMODELNAME = obj.get_val_raw(obj.meta_hash, 'ManufacturerModelName', true);
      [wip, statusWIP] = obj.get_val_raw(obj.meta_hash, 'sWipMemBlock', false);
      index = 0;

      switch lower(valSEQ)
        case {'b1v2d3d2', 'b1epi2d3d2'}
          index = 1;
        case {'b1epi4a3d2', 'b1epi2b3d2'}
          index = 2;
        case 'seste1d3d2'
          index = 14;
        case {'b1sev1a3d2' 'b1sev1b3d2' 'b1epi2f3d2' 'b1epi2g3d2' 'b1sev1a'}
          if contains(valMODELNAME,'Prisma','IgnoreCase',true)
            index = 14;
          elseif contains(valMODELNAME,'7T','IgnoreCase',true)
            index = 1;
          end
        otherwise
          warning(['%s (%s): Unknown sequence %s (%s). ' ...
                   'Unable to extract B1mapMixingTime'], ...
                  obj.filelist{obj.index}, fieldName, valSEQ, valPROT)
      end

      if statusWIP && index > 0
        if index <= numel(wip.alFree)
          value = wip.alFree(index) / 1000;
          result = true;
        else
          warning('%s (%s): Index %d out of range for sWipMemBlock', ...
                  obj.filelist{obj.index}, fieldName, index);
        end
      end
        
    otherwise
      [value, result] = obj.get_val_raw(obj.meta_hash, fieldName, true);
  end
end


function [res, status] = get_EchoSpacing(sequence, wip, model)
  res = [];
  status = false;

  switch lower(sequence)
    case {'b1v2d3d2', ... % 540 us - Prisma 
          'b1epi2b3d2', ... % 1mm protocol from WTCN'
          'b1epi2d3d2'}
      res = 0.540;
      status = true;
    case 'b1epi4a3d2' % 330 us - Allegra
      res = 0.330;
      status = true;
    case 'seste1d3d2'
      try
        res = (wip.alFree(5) * 2 + wip.alFree(6)) / 1000;
        status = true;
      catch
        warning('EchoSpacing: Unable to get WipParameters');
      end
   case {'b1sev1a3d2', 'b1sev1b3d2', 'b1epi2f3d2', 'b1epi2g3d2', 'b1sev1a'}
      % 7T and Prisma versions by Kerrin Pine
      if contains(model, '7T', 'IgnoreCase', true)
        res = 0.540;
        status = true;
      elseif contains(model, 'Prisma', 'IgnoreCase', true)
        try
          res = (wip.alFree(5) * 2 + wip.alFree(6)) / 1000;
          status = true;
        catch
          warning('EchoSpacing: Unable to get WipParameters');
        end
      end
   end
end

function [res, status] = get_EchoSpacing_BPP(PixelBandwidth)
  res = [];
  status = false;
  switch PixelBandwidth
    case 2300
      res = 0.540;
      status = true;
    case {3600, 3550}
      res = 0.330;
      status = true;
  end
end
