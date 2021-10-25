function [result, status] = get_val(obj, field, default_value)
  % retrieves value from currrent header
  % if value is not found, and default value is defined, 
  % returns default value

  result = [];
  status = false;

  switch obj.mode
    case 1
      [result, status] = obj.get_val_classic(field);
    case 2
      [result, status] = obj.get_val_classic(field);
    case 3
      [result, status] = obj.get_val_bids(field);
    otherwise
      warning('Undefined mode %d', obj.mode);
      result = [];
      status = false;
  end

  if ~status && exist('default_value', 'var')
    status = true;
    result = default_value;
  end
end
