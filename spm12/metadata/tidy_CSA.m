function tdyhdr = tidy_CSA(csahdr)
%==========================================================================
% function tdycsa = tidy_CSA(csahdr)
% DESCRIPTION:
% To rearrange CSAImageHeaderInfo, CSASeriesHeaderInfo, ... structures
% and make it a simplified, easier to browse structure (Siemens specific).
% USAGE:
% tdycsa = tidy_CSA(csahdr)
% where csahdr is a structure, content of the CSA field.
% WRITTEN BY: Evelyne Balteau - Cyclotron Research Centre
%==========================================================================

tdyhdr = [];
for ccsa = 1:length(csahdr)
    val = get_numaris4_val(csahdr,csahdr(ccsa).name);
    % if val is empty, let's not waste disk space with empty fields
    if ~isempty(val)
        % for elements having a value representation corresponding to a
        % numerical value (or an array of numerical values), convert char
        % into numbers. Here is a list of VR corresponding to numerical
        % values:  
        % - DS (Decimal String)
        % - FL (Floating Point Single)
        % - FD (Floating Point Double)
        % - IS (Integer String)
        % - OD (Other Double String)
        % - OF (Other Float String)
        % - SL (Signed Long)
        % - SS (Signed Short)
        % - UL (Unsigned Long)
        % - US (Unsigned Short)
        switch deblank(csahdr(ccsa).vr)
            case {'DS','FL','FD','IS','OD','OF','SL','SS','UL','US'}
                try
                    tmp = zeros(size(val,1),1);
                    for k = 1:size(val,1)
                        tmp(k)=str2num(val(k,:)); %#ok<ST2NM>
                    end
                    val = tmp;
                catch
                    fprintf(1,'Trouble reading CSA header %s (%s) = %s\n',csahdr(ccsa).name, deblank(csahdr(ccsa).vr), val(1,:));
                    val = [];
                end
            otherwise
                val = deblank(val);
        end
        % make sure no "-" in field name (invalid otherwise)
        csahdr(ccsa).name = strrep(csahdr(ccsa).name,'-','');
        tdyhdr.(csahdr(ccsa).name) = val;    
    end
end

%==========================================================================
% function val = get_numaris4_val(str,name)
% copied from spm_dicom_convert
%==========================================================================
function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str)
    if strcmp(deblank(str(i).name),name)
        if isfield(str(i),'nitems')
            for j=1:str(i).nitems
                if  str(i).item(j).xx(1)
                    val = [val {str(i).item(j).val}]; %#ok<*AGROW>
                end
            end
%         else
%             fprintf(1,'Found empty CSA header with name %s\n',name);
        end
        break;
    end
end
val = strvcat(val{:}); %#ok<DSTRVCT>