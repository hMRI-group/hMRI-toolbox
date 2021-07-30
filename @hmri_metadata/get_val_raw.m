function [res, status] = get_val_raw(inStruct, fieldName, casesensitive)
  % Recusively retrives value of the first found field in inStruct
  % If casesensitive is false, than field names are maches
  % are not case sensitive
  %
  % Arguments:
  % ----------
  %   inStruct: struct
  %       structure where the value of given name is searched for
  %   fieldName: char
  %       char array containing the field name to search, the leading
  %       and ending spaces are ignored
  %   casesensitive: bool
  %       if true, the field names are matched concidering case
  %
  % Returns:
  % --------
  %   res:
  %       found vaue, if corresponding field exists in structure,
  %       empty array [] otherwise
  %   status: bool
  %       true is searched field is found,
  %       false otherwise


  res = [];
  status = false;

  if ~isstruct(inStruct)
    return;
  end

  s1 = strtrim(fieldName);
  if ~casesensitive
    s1 = lower(s1);
  end

  f = fieldnames(inStruct(1));
  for i = 1:numel(f)
    if status
      break;
    end

    s2 = strtrim(f{i});
    if ~casesensitive
      s2 = lower(s2);
    end

    status = strcmp(s1,s2);

    if status
      res = inStruct(1).(f{i});
      break;
    end

    if iscell(inStruct(1).(f{i}))
      for cc = 1:numel(inStruct(1).(f{i}))
        [res, status] = hmri_metadata.get_val_raw(inStruct(1).(f{i}){cc}, ...
                                                  s1, casesensitive);
        if status
          break;
        end
      end
    end

    if isstruct(inStruct(1).(f{i}))
      [res, status] = hmri_metadata.get_val_raw(inStruct(1).(f{i}), ...
                                                s1, casesensitive);
      if status
        break;
      end
    end
  end

  if ischar(res)
    res = strtrim(res);
  end
end
