function res = find_bids_entity(fname, entity)
  % Split bids-formatted filename into parts and retrives
  % requested value from entity.

  [pth, base, ext] = fileparts(fname);
  base = [base, ext];
  res = [];

  l_entity = size(entity, 2);
  pos = 0;

  for i = 1:size(base,2)
    if base(i) == '_'
      sub = base(pos + 1: i - 1);
      pos = i;
      if strncmp(sub, entity, l_entity)
        res = sub(l_entity + 1: end);
        break
      end
    elseif base(i) == '.'
      if strcmp('suffix', entity)
        res = base(pos + 1: i - 1);
        break;
      elseif strcmp('extension', entity)
        res = base(i:end);
        break;
      end
    end
  end

end
