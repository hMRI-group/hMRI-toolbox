function [value, result] = get_val_bids(obj, fieldName)
  value = [];
  result = false;

  switch fieldName
    case 'RepetitionTime'
      % Valid for all vendors
      [value, result] = obj.get_val_raw(obj.meta_hash, 'RepetitionTime', true);
      if ~result
        [value, result] = obj.get_val_raw(obj.meta_hash, 'RepetitionTimeExcitation', true);
      end

      if result
        value = value * 1000; % converting to ms
      end

    case 'EchoTime'
      % Valid for all vendors
      [value, result] = obj.get_val_raw(obj.meta_hash, 'EchoTime', true);

      if result
        value = value * 1000; % converting to ms
      end

    case 'MT'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'MTstate', true);
      if result
        if ischar(value)
          if strcmpi(value, 'on')
            value = 1.;
          elseif strcmpi(value, 'off')
            value = 0.;
          else
            warning(['Invalid MTState value ' value]);
            value = 0.;
          end
        else
          value = double(value > 0);
        end
      else
        mt = find_bids_entity(obj.filelist{index}, 'mt-');
        if strcmpi(mt, 'on') || strcmp(mt, '1')
          result = true;
          value = 1.;
        elseif strcmpi(mt, 'off') || strcmp(mt, '0')
          result = true;
          value = 0.;
        end
      end

    case 'FieldStrength'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'MagneticFieldStrength', true);

    case 'PhaseEncodingDirectionSign'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'PhaseEncodingDirection', true);
      if result
        if value(end) == '-'
          value = -1;
        else
          value = 1;
        end
      else
        warning('PhaseEncodingDirection not found, using PhaseEncodingDirectionSign');
        [value, result] = obj.get_val_raw(obj.meta_hash, 'PhaseEncodingDirectionSign', true);
      end

    case 'PhaseEncodingDirection' % 'COL' (A>>P/P>>A/j) or 'ROW' (R>>L/L>>R/i)
      [value, result] = obj.get_val_raw(obj.meta_hash, 'PhaseEncodingDirection', true);
      if result
        if value == 'j'
          value = 'COL';
        elseif value == 'i'
          value = 'ROW'
        else
          warning(['Invalid PhaseEncodingDirection ' value]);
          value = [];
          result = false;
        end
      end

    case 'MultiBandFactor'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'MultibandAccelerationFactor', true);

    case 'isDWI'
        suffix = find_bids_entity(obj.filelist{obj.index}, 'suffix');
        result = true;
        if strcmpi(suffix, 'dwi')
          value = 1;
        else
          value = 0;
        end

    case 'epiReadoutDuration'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'epiReadoutTime', true);
      if ~result
        warning('epiReadoutTime not found using TotalReadoutTime');
        [value, result] = obj.get_val_raw(obj.meta_hash, 'TotalReadoutTime', true);
      end

      if result
        value = value * 1000;
      end

    case 'EchoSpacing'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'EffectiveEchoSpacing', true);
      if result
        value = value * 1000;
      end

    case 'B1mapMixingTime'
      [value, result] = obj.get_val_raw(obj.meta_hash, 'MixingTime', true);
      if result
        value = value * 1000;
      end

    otherwise
      [value, result] = obj.get_val_raw(obj.meta_hash, fieldName, true);
  end

  return;
end
