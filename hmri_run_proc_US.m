function out = hmri_run_proc_US(job)
% Deal with the spatial preprocessing, 1 subject at a time: segmentation of
% the MT and T1 images
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips
% but largely inspired by the batch from the past VBQ toolbox.

% Turning data organization from "per image type" into "per subject"
% because data are processed subject per subject.
% This relies alos on the previous toolbox, which allowed explicitly for a
% "per subject" setting of the data, so here we keep about the same code.
job = preproc_perimage_to_persubject(job);

% Initiliazign the output structure 'out'
% .tiss : struct-array with subfields
%           .c and .rc, for the native and Dartel imported
%           .wc and .mwc, for the warped and modulated
%          tissue class images
% .maps : struct-array with subfields 'wvols_pm' for the warped parametric
%         maps
% .def  : cell-array with the deformations for each subject.
for i=1:numel(job.tissue)
    out.tiss(i).c = {};
    out.tiss(i).rc = {};
    out.tiss(i).wc = {};
    out.tiss(i).mwc = {};
end
if numel(job.subjc(1).maps.vols_pm)
    for i=1:numel(job.subjc(1).maps.vols_pm)
        out.maps(i).wvols_pm = {};
    end
else
    out.maps.wvols_pm = {};
end
out.def.fn = {};

% looping over all the subjects.
for nm = 1:length(job.subjc)
    % Figure out where results are written out -> dn_output
    % and create it of needs be (per-subject option)
    same_dir = false;
    % pathes to struct image and parametric maps, could be different ones.
    struc_path = spm_file(job.subjc(nm).channel(1).vols{1},'path');
    if ~isempty(job.subjc(nm).maps.vols_pm)
        data_path = spm_file(job.subjc(nm).maps.vols_pm{1},'path');
    else
        data_path = struc_path;
    end
    if isfield(job.subjc(nm).output,'indir') && ...
            job.subjc(nm).output.indir == 1
        same_dir = true;
        dn_output = data_path;
    elseif isfield(job.subjc(nm).output,'outdir')
        dn_output = job.subjc(nm).output.outdir{1};
    elseif isfield(job.subjc(nm).output,'outdir_ps')
        % Get the subjects directory name, from data_ or struct_path???
        dn_subj = get_subject_dn(data_path); 
        dn_output = fullfile(job.subjc(nm).output.outdir_ps{1},dn_subj);
        if ~exist(dn_output,'dir')
            % Create subject sub-directory if necessary
            mkdir(dn_output);
        end
    end
    
    % Prepare and run 'spm_preproc' -> get tissue maps + deformation
    defsa.channel  = job.subjc(nm).channel;
    defsa.tissue   = job.tissue;
    defsa.warp     = job.warp;
    defsa.warp.vox = mean(job.many_sdatas.vox);
    defsa.warp.bb  = job.many_sdatas.bb;
    out.subjc(nm)  = spm_preproc_run(defsa);
    
    % Move segmentation output (if requested) and update 'out' structure:
    % all *c*.nii images, deformation field (y_*.nii), parameters 
    %  (*_seg8.mat), and bias corrections (m*.nii & BiasField_*.nii)
    if ~same_dir
        l_filter = {'^c[\d].*\.nii$','^rc[\d].*\.nii$', ...
            '^wc[\d].*\.nii$','^mwc[\d].*\.nii$','^y_.*\.nii$', ...
            '^iy_.*\.nii$','^.*_seg8.mat$','^BiasField_.*\.nii$'};
        f2move = spm_select('FPList',struc_path,l_filter);
        for ii=1:size(f2move,1)
            movefile(deblank(f2move(ii,:)),dn_output);
        end
        % Carefull with m*.nii files, which could be original
        % -> check channels
        if ~isempty(out.subjc(nm).channel)
            for ii = 1:numel(out.subjc(nm).channel)
                if ~isempty(out.subjc(nm).channel(ii).biascorr)
                    movefile(out.subjc(nm).channel(ii).biascorr{1},dn_output);
                end
            end
        end
        
        % Get the filenames updated
        out.subjc(nm) = update_path(out.subjc(nm),dn_output);
    end
    
    % Apply deformation on maps + get deformation map name
    defs.comp{1}.def = out.subjc(nm).fordef;
    % defs.ofname = '';
    if numel(job.subjc(nm).maps.vols_pm)
        defs.out{1}.pull.fnames = cellstr(char(char(job.subjc(nm).maps.vols_pm{:})));
        defs.out{1}.pull.savedir.saveusr{1} = dn_output;
        defs.out{1}.pull.interp = 1;
        defs.out{1}.pull.mask = 1;
        defs.out{1}.pull.fwhm = [0 0 0]; % no smoothing requester,
        % though at least vx_size/4 smoothing will still be applied!
        outdef = spm_deformations(defs);
    else
        outdef.warped = {};
    end
    
    % Save filenames as apropriate for 'out', keeping track of moved files!
    for i=1:numel(out.subjc(1).tiss)
        if isfield(out.subjc(nm).tiss(i), 'c')
            out.tiss(i).c = [out.tiss(i).c; out.subjc(nm).tiss(i).c];
        end
        if isfield(out.subjc(nm).tiss(i), 'rc')
            out.tiss(i).rc = [out.tiss(i).rc; out.subjc(nm).tiss(i).rc];
        end
        if isfield(out.subjc(nm).tiss(i), 'wc')
            out.tiss(i).wc = [out.tiss(i).wc; out.subjc(nm).tiss(i).wc];
        end
        if isfield(out.subjc(nm).tiss(i), 'mwc')
            out.tiss(i).mwc = [out.tiss(i).mwc; out.subjc(nm).tiss(i).mwc];
        end
    end
    for i=1:numel(outdef.warped)
        out.maps(i).wvols_pm{end+1,1} = outdef.warped{i};
    end
    out.def.fn{end+1,1} = defs.comp{1}.def{1};
end
end

% ========================================================================
%% SUBFUNCTIONS
% ========================================================================
function job = preproc_perimage_to_persubject(job)
% Rearrange data per subject for further preprocessing.

% Number of subjects from 1st set of structurals for segmentation
nSubj = numel(job.many_sdatas.channel(1).vols);
nChan = numel(job.many_sdatas.channel);
nPara = numel(job.many_sdatas.vols_pm);
for ii = 1:nSubj
    job.subjc(ii).output = job.many_sdatas.output;
    % Collect structurals
    job.subjc(ii).channel = job.many_sdatas.channel;
    for jj = 1:nChan
        job.subjc(ii).channel(jj).vols = job.many_sdatas.channel(jj).vols(ii);
    end
    % Collect parametric maps to warp, if any
    job.subjc(ii).maps.vols_pm = {};
    for kk = 1:nPara
        job.subjc(ii).maps.vols_pm{end+1,1} = ...
            job.many_sdatas.vols_pm{kk}{ii};
    end
end
end
%_______________________________________________________________________

function subjc_o = update_path(subjc_i,dn_output)
% Function to update the path of created files, when results are moved to
% another directory.

subjc_o = subjc_i; % At worst keep the same...

% Channel
if ~isempty(subjc_i.channel)
    % deal with 'biasfield' and 'biascorr' path
    for ii=1:numel(subjc_i.channel)
        subjc_o.channel(ii) = update_struct_path(subjc_i.channel(ii),dn_output);
    end
end

% Tissue
if ~isempty(subjc_i.tiss)
    for ii=1:numel(subjc_i.tiss)
        subjc_o.tiss(ii) = update_struct_path(subjc_i.tiss(ii),dn_output);
    end
end
    
% Parameters
if ~isempty(subjc_i.param)
    subjc_o.param = spm_file(subjc_i.param,'path',dn_output);
end

% Inverse deformation
if ~isempty(subjc_i.invdef)
    subjc_o.invdef = spm_file(subjc_i.invdef,'path',dn_output);
end

% Forward deformation
if ~isempty(subjc_i.fordef)
    subjc_o.fordef = spm_file(subjc_i.fordef,'path',dn_output);
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
