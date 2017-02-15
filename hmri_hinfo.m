function p = hmri_hinfo(P)
% Function to extract the TR/TE/FA values out of the descriptor field of 
% the nifti header

N = nifti(P);
nN = numel(N);
p(nN) = struct('tr',[],'te',[],'fa',[]);

hdr = get_metadata(P);
if ~isempty(hdr{1})
    for ii = 1:numel(N)
        p(ii).tr = get_metadata_val(hdr{ii},'RepetitionTime');
        p(ii).te = get_metadata_val(hdr{ii},'EchoTime')';
        p(ii).fa = get_metadata_val(hdr{ii},'FlipAngle');
    end
else
    for ii = 1:numel(N)
        tmp = regexp(N(ii).descrip,...
            'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
            'names');
        if isempty(tmp)
            p = [];
%             switch ii
%                 case 1
%                    prompt = strcat('Please enter the Echo Time (TE) of image #',...
%                         int2str(ii),' (',P(1,:),'): ');
%                    p(ii).te = input(prompt);
%                    prompt = strcat('Please enter the Repetition Time (TR) of image #',...
%                         int2str(ii),' (',P(1,:),'): ');
%                    p(ii).tr = input(prompt);
%                    prompt = strcat('Please enter the Flip Angle (FA) of image #',...
%                         int2str(ii),' (',P(1,:),'): ');
%                    p(ii).fa = input(prompt);
%                         otherwise
%                    prompt = strcat('Please enter the Echo Time (TE) of image #',...
%                         int2str(ii),' (',P(ii,:),'): ');
%                    p(ii).te = input(prompt);
%                     prompt = strcat('Please enter the Repetition Time (TR) of image #',...
%                         int2str(ii),' (',P(ii,:),') or use the same as before [TR = ',...
%                         num2str(p(ii-1).tr),']: ');
%                    p(ii).tr = input(prompt);
%                    if isempty(p(ii).tr)
%                      p(ii).tr = p(ii-1).tr;
%                    end
%                    prompt = strcat('Please enter the Flip Angle (FA) of image #',...
%                         int2str(ii),' (',P(ii,:),') or use the same as before [FA = ',...
%                         num2str(p(ii-1).fa),']: ');
%                    p(ii).fa = input(prompt);
%                    if isempty(p(ii).fa)
%                      p(ii).fa = p(ii-1).fa;
%                    end
%             end
            else
                p(ii).tr = str2num(tmp.tr);
                p(ii).te = str2num(tmp.te);
                p(ii).fa = str2num(tmp.fa);
        end
    end
end