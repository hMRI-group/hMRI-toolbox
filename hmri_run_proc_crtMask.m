function out = hmri_run_proc_crtMask(job)
% Function to run the mask creation routine
%
% The out structure has 2 fields:
% .fn_maskTC : filenames (char array) of the tissue specific masks
% .fn_meanTC : filenames (char array) of the mean tissue class images
%_______________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Centre

% Written by Christophe Phillips

% Unpack job input & sort out output directory
% n_TCs = numel(job.vols_smwc);     % #tissue classes
% for ii = 1:n_TCs
%     fn_smwTC{ii,1} = job.vols_smwc{ii};
% end
fn_smwTC = job.vols_smwc;
% options.val     = {threshTC noOverlap output};
opts = struct(...
    'minTCp', job.options.threshTC, ...
    'noOvl', job.options.noOverlap, ... 
    'outPth', '');
if isfield(job.options.output,'outdir') % use provided output dir
    outPth = job.options.output.outdir{1};
    if ~exist(outPth,'dir')
        mkdir(outPth); 
        fprintf('\nWARNING:\n');
        fprintf('\tCreating directory: %s\n\n',outPth)
    end
    opts.outPth = outPth;
end

% call function
[fn_maskTC,fn_meanTC] = hmri_proc_crtMask(fn_smwTC, opts);

% pack up output mask/mean image filenames
out.fn_maskTC = fn_maskTC;
out.fn_meanTC = fn_meanTC;

end

