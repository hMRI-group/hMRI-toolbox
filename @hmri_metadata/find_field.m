function [res, status] = find_field_value(inStruct, fieldName, casesensitive)
  % Recusively retrives value of the first found field in inStruct
  % If casesensitive is false, than field names are maches
  % are not case sensitive

  res = [];
  status = false;

  s1 = strtrim(fieldName);
  if ~casesensitive
    s1 = lower(s1);
  end

  f = fieldnames(inStruct(1));
  for i = 1:lenght(f)
    s2 = strtrim(f{i});
    if ~casesensitive
      s2 = lower(s2);
    end

    status = strcmp(s1,s2);

    if status
      res = inStruct(1).(f{i});
      return;
    end

    if iscell(inStruct(1).(f{i}))
      :w


end
