function tbx_wcomb_out = tbx_scfg_hmri_wcomb

%==========================================================================
% PURPOSE
% Calculation of robust combination of hMRI maps using error maps.
% This approach requires the acquisition of two successive datasets.
%==========================================================================

%% define variable: MPMs first run
in_vols1         = cfg_files;
in_vols1.tag     = 'in_vols1';
in_vols1.name    = 'MTsat, PD, and R1 maps from first run';
in_vols1.help    = {'Select MTsat, PD, and R1 maps from first run. The maps must be in the same order for the second run.'};
in_vols1.filter = 'image';
in_vols1.ufilter = '.*';
in_vols1.num     = [0 3];

%% define variable: MPMs second run
in_vols2         = cfg_files;
in_vols2.tag     = 'in_vols2';
in_vols2.name    = 'MTsat, PD, and R1 maps from second run';
in_vols2.help    = {'Select MTsat, PD, and R1 maps from second run. The maps must be in the same order for the first run.'};
in_vols2.filter = 'image';
in_vols2.ufilter = '.*';
in_vols2.num     = [0 3];

%% define variable: weights first run
in_weights1         = cfg_files;
in_weights1.tag     = 'in_weights1';
in_weights1.name    = 'Weight images for first run';
in_weights1.help    = {'Select weight images from first run corresponding to each map in same order as above.'};
in_weights1.filter = 'image';
in_weights1.ufilter = '.*param_error';
in_weights1.num     = [0 3];

%% define variable: weights second run
in_weights2         = cfg_files;
in_weights2.tag     = 'in_weights2';
in_weights2.name    = 'Weight images for second run';
in_weights2.help    = {'Select weight images from second run corresponding to each map in same order as above.'};
in_weights2.filter = 'image';
in_weights2.ufilter = '.*param_error';
in_weights2.num     = [0 3];

%% define variable: reference image
in_ref         = cfg_files;
in_ref.tag     = 'in_ref';
in_ref.name    = 'Reference image (or done for none)';
in_ref.help    = {'Select a reference image, to which all other data will be resampled. If no reference image is selected, the map from the first run will be used as reference.'};
in_ref.filter = 'image';
in_ref.ufilter = '.*';
in_ref.num     = [0 1];
in_ref.val     = {''};

%% define variable: mask image
in_msk         = cfg_files;
in_msk.tag     = 'in_msk';
in_msk.name    = 'Brain mask';
in_msk.help    = {'Select a brain mask image. This will improve the specifity of the robust combination approach.'};
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
                    'Robust combination of two MPMs from successive runs using error maps. This approach requires the acquisition of two successive runs of the MPM protocol. The proposed method and an example protocol is described in Mohammadi et al. 2022, NeuroImage'
};
tbx_wcomb_out.prog = @local_hmri_wcomb;
tbx_wcomb_out.vout = @out_hmri_wcomb;

end

%% dependencies
function out = local_hmri_wcomb(job)
nInputs = size(job.in_vols1,1);
Pout = cell(nInputs,1);
for n=1:nInputs
    Pout{n} = hmri_wcomb_2mpms(char(job.in_vols1(n)), char(job.in_vols2(n)), char(job.in_weights1(n)), char(job.in_weights2(n)), char(job.in_ref), char(job.in_msk));
end

extentwa = '_wa';
out.wafiles = spm_file(Pout, 'suffix', extentwa);

end

function dep = out_hmri_wcomb(job)
kk = 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Robust combination';
dep(kk).src_output = substruct('.','wafiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end