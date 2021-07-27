classdef hmri_metadata < handle
  properties
    filelist
    index
    mode
    meta_hash
  end

  methods
    function obj = hmri_metadata(files, varargin)
      % Basic object for retrieving metadata from MRI images
      % Supported formats:
      %     - nifti images with extended header
      %     - nifti images with json header
      %     - bids-formatted nifti images
      %
      % Arguments:
      % -----------
      %   files: char or cellstr
      %     file(s) to extract header
      %     
      % Parameters:
      % -----------
      %   'mode': char, one of --
      %       'automatic': metadata format will be retrieved automatically
      %       'classic': the hMRI header structure is assumed
      %       'header': extended nifti header assumed
      %       'bids': bid structure is assumed
      %

      if ischar(files)
        files = cellstr(files);
      end

      args = inputParser();
      args.addParameter('mode', 'automatic');
      args.parse(varargin{:});

      % Checking for existance of files
      for i = 1:numel(files)
        % remove ',1' from filename if needed
        res = regexp(files{i}, ',\d+$');
        if ~isempty(res)
          files{i} = files{i}(1:res - 1);
        end
        if ~exist(files{i}, 'file')
          error('File (%d) %s: not found', i, files{i});
        end
        [~, ~, ext] = fileparts(files{i});
        if ~strcmpi(ext, '.nii') and ~strcmpi(ext, '.json')
          error('File (%d) %s: unsupported extention %s', ...
                i, files{i}, ext);
        end
      end
      obj.filelist = files;

      switch args.Results.mode
        case 'automatic'
          obj.mode = 0;
        case 'classic'
          obj.mode = 1;
        case 'header'
          obj.mode = 2;
        case 'bids'
          obj.mode = 3;
        otherwise
          error('Unknown mode: %s', args.Results.mode);
      end

      % Loading first file
      obj.set_file(1, true);
    end

  end
end
