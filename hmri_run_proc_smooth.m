function out = hmri_run_proc_smooth(job)
% Function to run the smoothing/weighted averaging over a bunch of
% subjects, as defined in the batch interface.
% Data are selected in a 'many subject' style, i.e. all the images of one
% type are selected from many subjects at once!
% 
% The 'out' structure is organized as a structure out.tc where
% - tc is a cell-array of size {n_TCs x n_pams}
% - each element tc{ii,jj} is a cell array {n_subj x 1} with each subject's
%   smoothed data for the ii^th TC and jj^th MPM
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

% grab a few numbers, assuming data have been checked before
n_pams = numel(job.vols_pm);     % #parametric image types
n_TCs = numel(job.vols_tc);      % #tissue classes
n_subj = numel(job.vols_pm{1});  % #subjects

% disp([n_pams n_TCs n_subj])

% Find the list of tissue classes considered
l_TC = zeros(1,n_TCs);
for ii = 1:n_TCs
    tmp = regexp(job.vols_tc{ii}{1},'mwc(\d)', 'tokens');
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

for i_subj = 1:n_subj
    fn_wMPM = cell(n_pams,1);
    for jj = 1:n_pams
        fn_wMPM{jj} = job.vols_pm{jj}{i_subj};
    end
    fn_mwTC = cell(n_TCs,1);
    for jj = 1:n_TCs
        fn_mwTC{jj} = job.vols_tc{jj}{i_subj};
    end
    
    fn_out = hmri_proc_MPMsmooth(char(fn_wMPM), char(fn_mwTC), fn_TPM, job.fwhm, l_TC);
    
    for jj = 1:n_TCs
        for kk = 1:n_pams
            out.tc{jj,kk}{i_subj,1} = fn_out{kk}(jj,:); %#ok<*STRNU>
        end
    end
    
end

end

