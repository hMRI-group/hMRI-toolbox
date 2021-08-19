function [result, status] = get_list(obj, field, varargin)
  args = inputParser();
  args.addParameter('backup', '', @ischar);
  args.addParameter('default', [], @isnumeric);
  args.parse(varargin{:});

  result = [];
  status = false;

  [result, status] = obj.get_val(field);

  if ~status && ~isempty(args.Results.backup)
    old_index = obj.index;
    old_hash = obj.meta_hash;

    status = true;
    for i = 1:size(obj.filelist, 1)
      obj.set_file(i);
      [res, loc_status] = obj.get_val(args.Results.backup);

      if loc_status
        result(i) = res;
      else
        status = false;
      end
    end

    obj.index = old_index;
    obj.meta_hash = old_hash;
  end

  if ~isempty(args.Results.default)
    status = true;
  end

  if status && size(result, 1) ~= size(obj.filelist, 1)
    warning('%s: Retrieved value size mismatch number of files', field);
    status = false;
  end

end
