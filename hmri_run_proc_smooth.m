function out = hmri_run_proc_smooth(job)
% Function to run the smoothing/weighted averaging over a bunch of
% subjects, as defined in the batch interface.
% Data are selected in a 'many subject' style, i.e. all the images of one
% type are selected from many subjects at once!
% 
% The 'out' structure is organized as a structure with 2 fields
% .tc   : cell-array of size {n_TCs x n_pams}. Each element tc{ii,jj} is a 
%         cell array {n_subj x 1} with each subject's smoothed data for
%         the ii^th TC and jj^th MPM
% .smwc : cell-array of size {n_TCs x1}. Each element smwc{ii} is a char
%         array (n_subj x 1) with each subject's smooth modulated warped
%         ii^th tissue class
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

% grab a few numbers, assuming data have been checked before
n_pams = numel(job.vols_pm);     % #parametric image types
n_TCs = numel(job.vols_mwc);     % #tissue classes
n_subj = numel(job.vols_pm{1});  % #subjects
% disp([n_pams n_TCs n_subj])

% Define output folder for processing function
output = struct('outDir', [], 'option', 'same'); % Default -> keep same as input data
if isfield(job.output,'outdir') % -> everything in the same
    output.outDir = job.output.outdir{1};
    output.option = 'allin';
elseif isfield(job.output,'outdir_ps') % -> per suject organization
    output.outDir = job.output.outdir_ps{1};
    output.option = 'subjspec';
end

% Find the list of tissue classes considered
l_TC = zeros(1,n_TCs);
for ii = 1:n_TCs
    tmp = regexp(job.vols_mwc{ii}{1},'mwc(\d)', 'tokens');
    l_TC(ii) = str2num(tmp{1}{1}); %#ok<*ST2NM>
end 

% Get the TPM file, without any number
fn_TPM_i = cell(n_TCs,1);
for ii=1:n_TCs
    fn_TPM_i{ii} = spm_file(job.tpm{1},'number',l_TC(ii));
end
fn_TPM = char(fn_TPM_i);

% Loop over all the subjects and process them one at a time
out.tc = cell(n_TCs,n_pams);
out.smwc = cell(n_TCs,1);

for i_subj = 1:n_subj
    fn_wMPM = cell(n_pams,1);
    for jj = 1:n_pams
        fn_wMPM{jj} = job.vols_pm{jj}{i_subj};
    end
    [pth,~,~]   = fileparts(fn_wMPM{1});
    pth_out = get_output_path(pth,output);

    fn_mwTC = cell(n_TCs,1);
    for jj = 1:n_TCs
        fn_mwTC{jj} = job.vols_mwc{jj}{i_subj};
    end
    
    [fn_out, fn_smwc] = hmri_proc_MPMsmooth( ...
                char(fn_wMPM), char(fn_mwTC), fn_TPM, ...
                job.fwhm, l_TC, pth_out);
    
    for jj = 1:n_TCs
        for kk = 1:n_pams
            out.tc{jj,kk}{i_subj,1} = fn_out{kk}(jj,:); %#ok<*STRNU>
        end
        out.smwc{jj}{i_subj,1} = fn_smwc(jj,:);
    end
    
end

end
%__________________________________________________________________________

%__________________________________________________________________________
function pth = get_output_path(pth,output)

% Generate desired output path.
if ~isempty(output)
    switch output.option
        case 'same'
            % no change to output path
        case 'allin'
            % put everythin the same predefined folder
            pth = output.outDir;
        case 'subjspec'
            % keep per-subject organisation in predefined folder
            % and create it if necessary
            l_fsep = strfind(pth,filesep);
            lp_fsep = [0 l_fsep length(pth)+1];
            dn_subj = pth(lp_fsep(end-1)+1:lp_fsep(end)-1);
            pth = fullfile(output.outDir,dn_subj);
        otherwise
            % inconsistent specification -> no change to output path
            fprintf('\nWrong output path specification, use input data path.\n');
    end
    if ~exist(pth,'dir'), mkdir(pth); end
end

end
%__________________________________________________________________________

