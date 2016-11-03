function p = hinfo(P)
N = nifti(P);
for ii = 1:numel(N),
    tmp = regexp(N(ii).descrip,...
        'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
        'names');
    p(ii).tr = str2num(tmp.tr);
    p(ii).te = str2num(tmp.te);
    p(ii).fa = str2num(tmp.fa);
end
