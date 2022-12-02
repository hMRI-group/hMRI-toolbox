function tbx_wcomb_out = tbx_scfg_hmri_wcomb

%==========================================================================
% PURPOSE
% Calculation of weighted-combination of hMRI maps using error maps.
% This approach requires the acquisition of two successive datasets.
%==========================================================================

%% define variable: MPMs first run
in_vols1         = cfg_files;
in_vols1.tag     = 'in_vols1';
in_vols1.name    = 'MTsat, PD, and R1 maps for first run';
in_vols1.help    = {'Select MTsat, PD, and R1 maps for first run. The images should be in the same order for the second run.'};
in_vols1.filter = 'image';
in_vols1.ufilter = '.*';
in_vols1.num     = [0 3];

%% define variable: MPMs second run
in_vols2         = cfg_files;
in_vols2.tag     = 'in_vols2';
in_vols2.name    = 'MTsat, PD, and R1 maps for second run';
in_vols2.help    = {'Select MTsat, PD, and R1 maps for second run. The images should be in the same order for the first run.'};
in_vols2.filter = 'image';
in_vols2.ufilter = '.*';
in_vols2.num     = [0 3];

%% define variable: weights first run
in_weights1         = cfg_files;
in_weights1.tag     = 'in_weights1';
in_weights1.name    = 'Weight images for first run';
in_weights1.help    = {'Select three weight images for first run, corresponding to the three MPMs MT, PD, and T1.'};
in_weights1.filter = 'image';
in_weights1.ufilter = '.*param_error';
in_weights1.num     = [0 3];

%% define variable: weights second run
in_weights2         = cfg_files;
in_weights2.tag     = 'in_weights2';
in_weights2.name    = 'Weight images for second run';
in_weights2.help    = {'Select three weight images for second run, corresponding to the three MPMs MT, PD, and T1.'};
in_weights2.filter = 'image';
in_weights2.ufilter = '.*param_error';
in_weights2.num     = [0 3];
%% define variable: reference image
in_ref         = cfg_files;
in_ref.tag     = 'in_ref';
in_ref.name    = 'Refence image (or done for none)';
in_ref.help    = {'Select a reference image, to which all other data will be resampled to. In case no reference image is selected, the first MT image of the first run will be used as reference.'};
in_ref.filter = 'image';
in_ref.ufilter = '.*';
in_ref.num     = [0 1];
in_ref.val     = {''};
%% define variable: mask image
in_msk         = cfg_files;
in_msk.tag     = 'in_msk';
in_msk.name    = 'Brain mask';
in_msk.help    = {'Select a brain mask image. This will improve the specifity of the weighted-combination approach.'};
in_msk.filter = 'image';
in_msk.ufilter = '.*';
in_msk.num     = [0 1];
in_msk.val     = {''};
%% call local wcomb function
tbx_wcomb_out         = cfg_exbranch;
tbx_wcomb_out.tag     = 'tbx_scfg_hmri_wcomb';
tbx_wcomb_out.name    = 'Combine two successive MPM datasets';
tbx_wcomb_out.val     = {in_vols1 in_vols2 in_weights1 in_weights2 in_ref in_msk};
tbx_wcomb_out.help    = {
                    'Weighted-combination of two MPMs from successive runs using error maps. This approach requires the acquisition of two successive runs of the MPM protocol. The proposed method and an example protocol is described in Mohammadi et al. 2022, NeuroImage'
};
tbx_wcomb_out.prog = @local_hmri_wcomb;
tbx_wcomb_out.vout = @out_hmri_wcomb;
%% dependencies
function out = local_hmri_wcomb(job)
Pout = hmri_wcomb_2mpms(char(job.in_vols1), char(job.in_vols2), char(job.in_weights1), char(job.in_weights2), char(job.in_ref), char(job.in_msk));

extentwa = '_wa';
out.wafiles = spm_file(Pout, 'suffix', extentwa);

function dep = out_hmri_wcomb(job)
kk = 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Weighted Average';
dep(kk).src_output = substruct('.','wafiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
