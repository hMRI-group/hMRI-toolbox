function out = hmri_run_proc_US(job)
% Deal with the spatial preprocessing, 1 subject at a time: segmentation of
% the MT and T1 images
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips
% but largely inspired by the batch from the past VBQ toolbox.

% Turning data organization from "per image type" into "per subject"
% because data processed subject per subject.
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
for i=1:numel(job.subjc(1).maps.vols_pm)
    out.maps(i).wvols_pm = {};
end
out.def.fn = {};

% looping over all the subjects.
for nm = 1:length(job.subjc)
    % Prepare and run 'spm_preproc' -> get tissue maps + deformation
    defsa.channel = job.subjc(nm).struct(1);
    defsa.channel.vols = job.subjc(nm).struct(1).s_vols;
    defsa.tissue  = job.tissue;
    defsa.warp    = job.warp;
    out.subjc(nm) = spm_preproc_run(defsa);
    
    % Apply deformation on strcut/maps + build deformation map
    defs.comp{1}.def = spm_file(job.subjc(nm).struct(1).s_vols, ...
        'prefix', 'y_', 'ext', '.nii'); % def map fname
    % defs.ofname = '';
    defs.out{1}.pull.fnames = cellstr(char(char(job.subjc(nm).maps.vols_pm{:})));
    if isfield(job.subjc(nm).output,'indir') && job.subjc(nm).output.indir == 1
        defs.out{1}.pull.savedir.saveusr{1} = ...
            spm_file(job.subjc(nm).maps.vols_pm{1},'path');
    else
        defs.out{1}.pull.savedir.saveusr{1} = job.subjc(nm).output.outdir{1};
    end
    defs.out{1}.pull.interp = 1;
    defs.out{1}.pull.mask = 1;
    defs.out{1}.pull.fwhm = [0 0 0]; % no smoothing requester,
    % though at least vx_size/4 smoothing will still be applied!
    outdef = spm_deformations(defs);
    
    % Save filenames as apropriate for 'out'
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
for i = 1:numel(job.many_sdatas.rstruct.s_vols)
    job.subjc(i).output = job.many_sdatas.output;
    job.subjc(i).struct = job.many_sdatas.rstruct;
    job.subjc(i).struct.s_vols = ...
        job.many_sdatas.rstruct.s_vols(i);
    job.subjc(i).maps.vols_pm = {};
    for k = 1:numel(job.many_sdatas.vols_pm)
        job.subjc(i).maps.vols_pm{end+1,1} = ...
            job.many_sdatas.vols_pm{k}{i};
    end
end
end
%_______________________________________________________________________
