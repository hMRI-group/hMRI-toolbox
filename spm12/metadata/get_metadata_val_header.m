function [parValue, parLocation] = get_metadata_val_header(varargin)
  % there seems to be no metadata available (mstruc is not a structure
  % and is likely to be empty). Still try a few tricks in case either
  % TR/TE/FA is the target field and available in the description field:
  inParName = varargin{2};
  if ischar(varargin{1}) && strcmp(spm_file(varargin{1},'ext'),'nii')
    fnam = varargin{1};
    fprintf(1,'\nWARNING: No metadata available in %s.\n', fnam);
    N = nifti(fnam);
    p = struct('tr',[],'te',[],'fa',[]);
    tmp = regexp(N.descrip,'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg','names');
    if ~isempty(tmp)
      p.tr = str2double(tmp.tr);
      p.te = str2double(tmp.te);
      p.fa = str2double(tmp.fa);
      switch inParName
        case 'RepetitionTime' % [ms]
          cRes = 1;
          parLocation{cRes} = 'NiftiDescriptionField';
          parValue{cRes} = p.tr;
          fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f ms\n', inParName, parValue{1});
            
        case 'EchoTime'
          cRes = 1;
          parLocation{cRes} = 'NiftiDescriptionField';
          parValue{cRes} = p.te;
          fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f ms\n', inParName, parValue{1});
            
        case 'FlipAngle'
          cRes = 1;
          parLocation{cRes} = 'NiftiDescriptionField';
          parValue{cRes} = p.fa;
          fprintf(1,'Parameter available in NIFTI description field:\n\t%s = %5.2f deg\n', inParName, parValue{1});
            
        otherwise
          parLocation = [];
          parValue = [];
          fprintf(1,'Not able to retrieve any value for parameter %s.\n', inParName);
      end
    end
  else
      error('Invalid input arguments. Type help get_metadata_val for correct syntax.');
  end
end
