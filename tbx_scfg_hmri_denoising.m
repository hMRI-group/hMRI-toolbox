function denoise = tbx_scfg_hmri_denoising

%--------------------------------------------------------------------------
% To enable/disable pop-up messages for all warnings - recommended when
% piloting the data processing.
%--------------------------------------------------------------------------
popup        = cfg_menu;
popup.tag    = 'popup';
popup.name   = 'Pop-up warnings';
popup.help   = {['The user can review and keep track of all the information ' ...
    'collected, including warnings and other messages coming up during ' ...
    'the creation of the maps. By default, the information is logged in ' ...
    'the Matlab Command Window, in a log file saved in the "Results" ' ...
    'directory, and when more critical, displayed as a pop-up message.'], ...
    ['The latter must be disabled for processing series of datasets (since it ' ...
    'blocks the execution of the code) but it is strongly recommended to ' ...
    'leave it enabled when piloting the data processing (single subject) ' ...
    'to read through and acknowledge every message and make sure ' ...
    'everything is set up properly before running the processing on a ' ...
    'whole group.'], ...
    ['More information about the various messages and action to be taken ' ...
    '(or not) accordingly can be found on the hMRI-Toolbox WIKI (http://hmri.info). ' ...
    'In particular, see the "Debug tips & tricks" section.']};
popup.labels = {'Disable' 'Enable'};
popup.values = {false true};
popup.val = {true};

% ---------------------------------------------------------------------
% menu denoisingtype
% ---------------------------------------------------------------------
denoisingtype = tbx_scfg_hmri_denoising_menu;

% ---------------------------------------------------------------------
% indir Input directory as output directory
% ---------------------------------------------------------------------
indir         = cfg_entry;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {['Output files will be written to the same folder ' ...
    'as each corresponding input file.']};
indir.strtype = 's';
indir.num     = [1 Inf];
indir.val     = {'yes'};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where output files will be written to.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as the input ' ...
    'directory for each input file or user selected']};
output.values  = {indir outdir};
output.val = {indir};

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj            = cfg_branch;
subj.tag        = 'subj';
subj.name       = 'Subject';
subj.help       = {'Specify a subject for denoising.'};
subj.val        = {output denoisingtype popup};

% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
sdata           = cfg_repeat;
sdata.tag       = 'data';
sdata.name      = 'Few Subjects';
sdata.help      = {'Specify the number of subjects.'};
sdata.num       = [1 Inf];
sdata.val       = {subj};
sdata.values    = {subj};

% ---------------------------------------------------------------------
% denoise images
% ---------------------------------------------------------------------
denoise         = cfg_exbranch;
denoise.tag     = 'denoise';
denoise.name    = 'Denoising';
denoise.val     = { sdata };
denoise.help    = {'Denoising of raw/processed images with different methods'};
denoise.prog    = @hmri_run_denoising;
denoise.vout    = @vout_create;

end

% ========================================================================
%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_denoising(job)
% is defined separately.
%_______________________________________________________________________

function dep = vout_create(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated, depending on the denoising type.

dnfield = fieldnames(job.subj.denoisingtype);
denoisingmethod = dnfield{1};

switch denoisingmethod
    case 'lcpca_denoise'
        % define variables and initialize cfg_dep based on availibility of phase images
        arrayLength = numel(job.subj.denoisingtype.lcpca_denoise.mag_input);
        % phase_bool= any(~cellfun(@isempty, job.subj.denoisingtype.lcpca_denoise.phase_input));
        phase_bool = isempty(job.subj.denoisingtype.lcpca_denoise.phase_input);
        if phase_bool
            cdep(1,2*arrayLength) = cfg_dep;
        else
            cdep(1,arrayLength) = cfg_dep;
        end

        % iterate to generate dependency tags for outputs
        for i=1:numel(job.subj)
            for k =1:2*arrayLength
                if k<=arrayLength
                    cdep(k)            = cfg_dep;
                    cdep(k).sname      = sprintf('lcpcaDenoised_magnitude%d',k);
                    idxstr = ['DenoisedMagnitude' int2str(k)];
                    cdep(k).src_output = substruct('.','subj','()',{i},'.',idxstr,'()',{':'});
                    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

                elseif k>arrayLength && ~phase_bool
                    cdep(k)            = cfg_dep;
                    cdep(k).sname      = sprintf('lcpcaDenoised_phase%d',k-arrayLength);
                    idxstr = ['DenoisedPhase' int2str(k-arrayLength)];
                    cdep(k).src_output = substruct('.','subj','()',{i},'.',idxstr,'()',{':'});
                    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                else
                    break
                end
            end
        end
        dep = cdep;
end
end