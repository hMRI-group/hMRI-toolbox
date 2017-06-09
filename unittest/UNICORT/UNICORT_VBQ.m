% List of open inputs
% Multiparameter & UNICORT_B1 images: Output directory - cfg_files
nrun = 1; % enter the number of runs here
jobfile = {'/data/pt_phy048/SD/hMRI-example-data-master-with-UnitTest/Leipzig_dataset/Prisma/maps/actualTB/UNICORT_VBQ_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Multiparameter & UNICORT_B1 images: Output directory - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
