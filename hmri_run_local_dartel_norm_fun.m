function out = hmri_run_local_dartel_norm_fun(job)
% function out = hmri_dartel_norm_fun_local(job)
%  derived from spm_dartel_norm_fun_local

% TO DO:
%  Build the output 'out' structure with created file names! 

out = struct();

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
spm_dartel_norm_fun(feds);

% MFC - Now take the MPMs and do a regular normalisation but don't apply
% Jacobian modulation. Produces ws* images.
for mm=1:length(job.subjd)
    feds.data.subj(mm).images = job.subjd(mm).mp_vols;
end
feds.preserve = 0;
spm_dartel_norm_fun(feds);

% Specify the out structure with created file names

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

