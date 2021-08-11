function set_file(obj, index, force)
  % Set current file to given index
  % 
  % Arguments:
  % ----------
  %   index: int
  %     index of file to load, must be within allowed range
  %   force: boolean, optional
  %     if true, the file header is loaded even if header is 
  %     already loaded

  if ~exist('force', 'var')
    force = false;
  end
  fname = obj.filelist{index};
  [path, basename, ext] = fileparts(fname);
  if ~force && strcmp(fname, obj.filelist{obj.index})
    return;
  end

  obj.index = index;
  json_file = fullfile(path, [basename '.json']);
  if exist(json_file, 'file')
    has_json = true;
  else
    has_json = false;
  end

  switch obj.mode
    case 0 % automatic mode
      if has_json
        hdr = spm_jsonread(json_file);
        if isfield(hdr, 'acqpar')
          hdr = hdr.acqpar;
          size(hdr)
          obj.mode = 1;
        else
          obj.mode = 3;
        end
      else
        fid = fopen(fname,'r');
        % skeeping standard header
        fseek(fid, 348);
        isHdrExtended = fread(fid, 4, 'uint8');
        if isHdrExtended(1)
          ext_hdr_size = fread(fid, 1, 'uint32');
          % read the ID of the extension (32-bit number) = anything > 4,
          % arbitraryly set to 27 for now, just because I like this number
          % Comment of original coder
          ext_hdr_id = fread(fid, 1, 'uint32'); 
          jsonhdr = fscanf(fid, '%c', ext_hdr_size-8);
          hdr = spm_jsonread(jsonhdr);
          obj.mode = 2;
        else
          error('File (%d) %s: only extended header are supported', ...
                i, [basename, ext]);
        end
      end
    case 1 % classic hdr
      if has_json
        hdr = spm_jsonread(json_file);
        if isfield(hdr, 'acqpar')
          hdr = hdr.acqpar;
        else
          error('File (%d) %s: missing ''acqpar'' field', ...
                i, [basename, ext]);
        end
      else
        error('File (%d) %s: missing json file', ...
              i, [basename, ext]);
      end

    case 2
        fid = fopen(fname,'r');
        % skeeping standard header
        fseek(fid, 348);
        isHdrExtended = fread(fid, 4, 'uint8');
        if isHdrExtended(1)
          ext_hdr_size = fread(fid, 1, 'uint32');
          % read the ID of the extension (32-bit number) = anything > 4,
          % arbitraryly set to 27 for now, just because I like this number
          % Comment of original coder
          ext_hdr_id = fread(fid, 1, 'uint32'); 
          jsonhdr = fscanf(fid, '%c', ext_hdr_size-8);
          hdr = spm_jsonread(jsonhdr);
        else
          error('File (%d) %s: only extended header are supported', ...
                i, [basename, ext]);
        end

    case 3
      if has_json
        hdr = spm_jsonread(json_file);
        if isfield(hdr, 'acqpar')
          error('File (%d) %s: missing ''acqpar'' field', ...
                i, [basename, ext]);
        end
      else
        error('File (%d) %s: missing json file', ...
              i, [basename, ext]);
      end
  end
  obj.meta_hash = hdr;
  obj.index = index;
end
