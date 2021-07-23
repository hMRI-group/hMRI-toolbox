classdef hmri_metadata
  properties
    filelist
    index
    mode
    override
    use_defaults
    meta_hash
  end

  methods
    function obj = hmri_metadata(files, varargin)
      if ischar(files)
        files = cellstr(files);
      end

      args = inputParser();
      args.addParameter('mode', 'automatic');
      args.addParameter('ovverride', []);
      args.parse(varargin{:});

      % Checking for existance of files
      for i = 1:numel(files)
        % TODO remove ',1' from filename if needed
        res = regexp(files{i}, ',\d+$');
        if ~empty(res)
          files{i} = files{i}(1:res - 1);
        end
        if ~exist(files{i}, 'file')
          error('File (%d) %s: not found', i, files{i});
        end
        [~, ~, ext] = fileparts(files{i});
        if ~strcmpi(ext, '.nii') and ~~strcmpi(ext, '.json')
          error('File (%d) %s: unsupported extention %s', ...
                i, files{i}, ext);
        end
      end

      switch args.Results.mode
        case 'automatic'
          self.mode = 0
        case 'classic'
          self.mode = 1
        case 'header'
          self.mode = 2
        case 'bids'
          self.mode = 3
        otherwise
          error('Unknown mode: %s', args.Results.mode);
      end

      if isempty(args.Results.ovverride)
        obj.override = [];
      elseif ~isstruct(args.Results.ovverride)
        error('Ovverride parameter must be either struct or empty');
      else
        obj.override = args.Results.ovverride;
      end

      % Loading first file
      obj.set_file(1, true);
    end

    function set_file(index, force)
      if ~exist(force, 'var')
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
              isBids = false;
            else
              isBids = true;
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
              isBids = false;
            else
              error('File (%d) %s: only extended header are supported', ...
                    i, files{i});
            end
          end
        case 1 % classic hdr
          if has_json
            hdr = spm_jsonread(json_file);
            if isfield(hdr, 'acqpar')
              isBids = false;
            else
              error('File (%d) %s: missing ''acqpar'' field', ...
                    i, files{i});
            end
          else
            error('File (%d) %s: missing json file', ...
                  i, files{i});
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
              isBids = false;
            else
              error('File (%d) %s: only extended header are supported', ...
                    i, files{i});
            end

        case 3
          if has_json
            hdr = spm_jsonread(json_file);
            if isfield(hdr, 'acqpar')
              error('File (%d) %s: missing ''acqpar'' field', ...
                    i, files{i});
            else
              isBids = true;
            end
          else
            error('File (%d) %s: missing json file', ...
                  i, files{i});
          end

      end
      
      if isBids 
    end
end
