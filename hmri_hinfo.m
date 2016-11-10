function p = hmri_hinfo(P)
% Function to extract the TR/TE/FA values out of the descriptor field of 
% the nifti header

N = nifti(P);
nN = numel(N);
p(nN) = struct('tr',[],'te',[],'fa',[]);
for ii = 1:numel(N)
    tmp = regexp(N(ii).descrip,...
        'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
        'names');
    p(ii).tr = str2num(tmp.tr);
    p(ii).te = str2num(tmp.te); %#ok<*ST2NM>
    p(ii).fa = str2num(tmp.fa);
end

end