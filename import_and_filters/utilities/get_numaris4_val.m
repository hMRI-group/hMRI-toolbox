%==========================================================================
% function val = get_numaris4_val(str,name)
% from spm_dicom_convert (SPM12)
% $Id: spm_dicom_convert.m 6773 2016-04-20 09:26:59Z john $
%==========================================================================
function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str)
    if strcmp(deblank(str(i).name),name)
        for j=1:str(i).nitems
            if  str(i).item(j).xx(1)
                val = [val {deblank(str(i).item(j).val)}];
            end
        end
        break;
    end
end
val = strvcat(val{:});
