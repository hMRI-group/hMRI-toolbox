function out = hmri_run_proc_dartel_norm(job)
% function out = hmri_run_proc_dartel_norm(job)
% derived from spm_dartel_norm_fun_local
% 
% It applies the "Dartel - Normalize to MNI" onto the tissue classes and
% parametric maps, from the job created in the batch.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

%  Build the output 'out' structure with created filenames in cell arrays
out.vols_mwtc = {}; % filenames are organised as {1 x #TCs} of {#subj x 1}
out.vols_wpm = {};  % filenames are organised as {1 x #PMs} of {#subj x 1} 


% MFC - Setting up feds structure which has a copy of the (reordered) subject info:
job = dartel_perimage_to_persubject(job);

feds.template = job.template;
feds.vox      = job.vox;
feds.bb       = job.bb;
feds.fwhm     = [0 0 0];
for m = 1:length(job.subjd)
    feds.data.subj(m).flowfield = job.subjd(m).flowfield;
    feds.data.subj(m).images    = job.subjd(m).tc_vols;
end

% MFC - Jacobian modulation correction to preserve total signal intensity:
feds.preserve = 1;

% MFC - Produces mwc* images, i.e. modulated, spatially normalised images.
% This produces w = |Dphi|t(phi), the product of the Jacobian determinants
% of deformation phi and the tissue class image warped by phi, as per
% Draganski 2011, NeuroImage.
out_tc = spm_dartel_norm_fun(feds);

% Specify the output with created file names, based on out_tc
out.vols_mwtc  = cell(1,numel(job.subjd(1).tc_vols));
for ii=1:numel(job.subjd(1).tc_vols)
    out.vols_mwtc{ii}  = cell(numel(job.subjd),1);
    for jj=1:numel(job.subjd)
        out.vols_mwtc{ii}{jj} = out_tc{jj}{ii};
    end
end

% MFC - Now take the MPMs and do a regular normalisation but don't apply
% Jacobian modulation. Produces ws* images.
for mm = 1:length(job.subjd)
    feds.data.subj(mm).images = job.subjd(mm).mp_vols;
end
feds.preserve = 0;
out_pm = spm_dartel_norm_fun(feds);

% Specify the output with created file names, based on out_pm
out.vols_wpm  = cell(1,numel(job.subjd(1).mp_vols));
for ii = 1:numel(job.subjd(1).mp_vols)
    out.vols_wpm{ii} = cell(numel(job.subjd),1);
    for jj = 1:numel(job.subjd)
        out.vols_wpm{ii}{jj} = out_pm{jj}{ii};
    end
end

end

% ========================================================================
%% SUBFUNCTIONS
% ========================================================================

function job = dartel_perimage_to_persubject(job)
% Rearrange data per subject for dartel processing.
for i=1:numel(job.multsdata.vols_field)
    job.subjd(i).tc_vols = {};
    for j=1:numel(job.multsdata.vols_tc)
        job.subjd(i).tc_vols{j} = job.multsdata.vols_tc{j}{i};
    end
    job.subjd(i).flowfield = {};
    job.subjd(i).flowfield{1} = job.multsdata.vols_field{i};
    job.subjd(i).mp_vols = {};
    for j=1:numel(job.multsdata.vols_pm)
        job.subjd(i).mp_vols{j} = job.multsdata.vols_pm{j}{i};
    end
end
end
%_______________________________________________________________________

