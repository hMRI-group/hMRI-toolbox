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
out.vols_mwc = {}; % filenames are organised as {1 x #mwc} of {#subj x 1}
out.vols_wpm = {};  % filenames are organised as {1 x #pm} of {#subj x 1}

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

% Jacobian modulation correction to preserve total signal intensity:
feds.preserve = 1;

% Produce mwc* images, i.e. modulated, spatially normalised images.
% This produces w = |Dphi|t(phi), the product of the Jacobian determinants
% of deformation phi and the tissue class image warped by phi, as per
% Draganski 2011, NeuroImage.
out_tc = spm_dartel_norm_fun(feds);

% Specify the output with created file names, based on out_tc
out.vols_mwc  = cell(1,numel(job.subjd(1).tc_vols));
for ii=1:numel(job.subjd(1).tc_vols)
    out.vols_mwc{ii}  = cell(numel(job.subjd),1);
    for jj=1:numel(job.subjd)
        out.vols_mwc{ii}{jj} = out_tc{jj}{ii};
    end
end

% Now take the MPMs and do a regular normalisation but don't apply
% Jacobian modulation. Produces ws* images.
feds.preserve = 0;
for mm = 1:length(job.subjd)
    feds.data.subj(mm).images = job.subjd(mm).mp_vols;
end
out_pm = spm_dartel_norm_fun(feds);

% Move things and update path of saved filenames
% NOTE: to be removed if this ends up part of the native Dartel function!
if isfield(job.output,'indir') && ...
        job.output.indir == 1
    % leave as it is and do nothing. :-)
else
    for ii=1:numel(out_pm) % for each subject
        % find output directory
        if isfield(job.output,'outdir') % -> everything in the same
            dn_output = job.output.outdir{1};
        elseif isfield(job.output,'outdir_ps')
            % Get the subjects directory name, from maps
            dn_subj = get_subject_dn(job.subjd(1).mp_vols{1});
            dn_output = fullfile(job.output.outdir_ps{1},dn_subj);
            if ~exist(dn_output,'dir')
                % Create subject sub-directory if necessary
                mkdir(dn_output);
            end
        else
            warning('hmri:norm2mni','Incorrect output directory specification.');
            break
        end
        % move files
        for jj=1:numel(out_pm{ii})
            movefile(out_pm{ii}{jj},dn_output)
        end
        % update path
        out_pm{ii} = spm_file(out_pm{ii},'path',dn_output);
    end
end

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
        job.subjd(i).tc_vols{j,1} = job.multsdata.vols_tc{j}{i};
    end
    job.subjd(i).flowfield = {};
    job.subjd(i).flowfield{1} = job.multsdata.vols_field{i};
    job.subjd(i).mp_vols = {};
    for j=1:numel(job.multsdata.vols_pm)
        job.subjd(i).mp_vols{j,1} = job.multsdata.vols_pm{j}{i};
    end
end
end
%_______________________________________________________________________


function st_o = update_struct_path(st_i,dn_o)
% Small function to update the path of filenames stored in subfield of an
% input structure 'st_i'.

field_nm = fieldnames(st_i);
st_o = st_i;
for ii = 1:numel(field_nm)
    if ~isempty(st_i.(field_nm{ii}))
        st_o.(field_nm{ii}) = spm_file(st_i.(field_nm{ii}),'path',dn_o);
    end
end
end
%_______________________________________________________________________

function dn_subj = get_subject_dn(data_path)
% Extract the subject's directory name from the path to its data

% Fist split the path into its sub-parts
l_fsep = strfind(data_path,filesep);
n_fsep = numel(l_fsep);
lp_fsep = [0 l_fsep length(data_path)+1];
pth_parts = cell(1,n_fsep);
for ii=1:(n_fsep+1)
    pth_parts{ii} = data_path(lp_fsep(ii)+1:lp_fsep(ii+1)-1);
end

% Pick up last one
dn_subj = pth_parts{end};

end
%_______________________________________________________________________