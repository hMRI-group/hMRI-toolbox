function out = hmri_run_local_preproc(job)
% Deal with the spatial preprocessing, 1 subject at a time: segmentation of
% the MT and T1 images
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id$

job = preproc_perimage_to_persubject(job);

for i=1:numel(job.tissue)
    out.tiss(i).c = {};
    out.tiss(i).rc = {};
end

for i=1:numel(job.subjc(1).maps.mp_vols)
    out.maps(i).mp_vols = {};
end

for nm = 1:length(job.subjc)
    defsa.channel = job.subjc(nm).struct(1);
    defsa.channel.vols = job.subjc(nm).struct(1).s_vols;
    defsa.tissue  = job.tissue;
    defsa.warp    = job.warp;
    out.subjc(nm) = spm_preproc_run(defsa);
    defs.comp{1}.def = ...
        strcat(spm_str_manip(job.subjc(nm).struct(1).s_vols,'h'), ...
        filesep,'y_',spm_str_manip(job.subjc(nm).struct(1).s_vols,'tr'),'.nii');
    % defs.ofname = '';
    defs.out{1}.pull.fnames = cellstr(char(char(job.subjc(nm).maps.mp_vols{:})));
    if isfield(job.subjc(nm).output,'indir') && job.subjc(nm).output.indir == 1
        defs.out{1}.pull.savedir.saveusr{1} = ...
            spm_str_manip(job.subjc(nm).maps.mp_vols{1},'h');
    else
        defs.out{1}.pull.savedir.saveusr{1} = job.subjc(nm).output.outdir{1};
    end
    defs.out{1}.pull.interp = 1;
    defs.out{1}.pull.mask = 1;
    defs.out{1}.pull.fwhm = [0 0 0];
    outdef = spm_deformations(defs);
    
    for i=1:numel(out.subjc(1).tiss)
        if isfield(out.subjc(nm).tiss(i), 'c')
            out.tiss(i).c = [out.tiss(i).c; out.subjc(nm).tiss(i).c];
        end
        if isfield(out.subjc(nm).tiss(i), 'rc')
            out.tiss(i).rc = [out.tiss(i).rc; out.subjc(nm).tiss(i).rc];
        end
    end
    for i=1:numel(outdef.warped)
        out.maps(i).mp_vols{end+1} = outdef.warped{i};
    end
    
    % Create sum of GM/WM
    for i=1:length(outdef.warped)
        if isfield(job.subjc(nm).output,'indir') && job.subjc(nm).output.indir == 1
            p = spm_str_manip(job.subjc(nm).maps.mp_vols{1},'h'); %#ok<*NASGU>
        else
            p = job.subjc(nm).output.outdir{1};
        end
        c1 = insert_pref(job.subjc(nm).struct(1).s_vols{1},'mwc1'); % modulated warped GM
        c2 = insert_pref(job.subjc(nm).struct(1).s_vols{1},'mwc2'); % modulated warped WM
        f  = insert_pref(job.subjc(nm).maps.mp_vols{i},'w');  % removed s f=outdef.warped{i};
        c  = spm_imcalc(char(char(c1),char(c2)),insert_pref(f,'bb_'),'(i1+i2)');
        c  = c.fname; % file with sum of GM & WM
        m_c1 = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',1']; % GM tpm
        m_c2 = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',2']; % WM tpm
        m_c  = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',6']; % Outside head tpm
        p1 = spm_imcalc(char(char(c1),char(f),m_c1),insert_pref(f,'p1_'), ...
            '(i1.*i2).*(i3>0.05)');
        p1 = p1.fname; % MP weighted with its own GM, and a priori GM>.05
        p2 = spm_imcalc(char(char(c2),char(f),m_c2),insert_pref(f,'p2_'), ...
            '(i1.*i2).*(i3>0.05)');
        p2 = p2.fname; % MP weighted with its own WM, and a priori GM>.05
        pp = spm_imcalc(char(char(c),char(f),m_c),insert_pref(f,'p_'), ...
            '(i1.*i2).*((1-i3)>0.05)');
        pp = pp.fname; % MP weighted with its own GM+WM, and a priori in the head > .05
        m1 = insert_pref(c1,'s'); spm_smooth(c1,m1,job.fwhm); % smooth mwc1
        m2 = insert_pref(c2,'s'); spm_smooth(c2,m2,job.fwhm); % smooth mwc2
        m  = insert_pref(c,'s');  spm_smooth(c,m,job.fwhm); % smooth mwc1+mwc2
        n1 = insert_pref(p1,'s'); spm_smooth(p1,n1,job.fwhm); % smooth weighted(GM) MP
        n2 = insert_pref(p2,'s'); spm_smooth(p2,n2,job.fwhm); % smooth weighted(WM) MP
        n  = insert_pref(pp,'s'); spm_smooth(pp,n,job.fwhm); % smooth weighted(GM+WM) MP
        q1 = spm_imcalc(char(n1,m1,m1),insert_pref(p1,'fin_uni_'), ...
            '(i1./i2).*(i3>0.05)'); % Signal as in paper
        q1 = q1.fname;
        q2 = spm_imcalc(char(n2,m2,m2),insert_pref(p2,'fin_uni_'), ...
            '(i1./i2).*(i3>0.05)');
        q2 = q2.fname;
        q  = spm_imcalc(char(n,m,m),insert_pref(pp,'fin_uni_bb_'), ...
            '(i1./i2).*((i3)>0.05)');
        q  = q.fname;
        delfiles = strrep({c,p1,p2,m1,m2,n1,n2,pp,m,n,q},'.nii,1','.nii');
        for ii=1:numel(delfiles)
            delete(delfiles{ii});
        end
    end
    
end
end

% ========================================================================
%% SUBFUNCTIONS
% ========================================================================
function job = preproc_perimage_to_persubject(job)
% Rearrange data per subject for further preprocessing.
if isfield(job, 'many_few_sdatas')
    if isfield(job.many_few_sdatas, 'subjc')
        job.subjc = job.many_few_sdatas.subjc;
    else
        for i = 1:numel(job.many_few_sdatas.many_sdatas.struct.s_vols)
            job.subjc(i).output = job.many_few_sdatas.many_sdatas.output;
            job.subjc(i).struct = job.many_few_sdatas.many_sdatas.struct;
            job.subjc(i).struct.s_vols = ...
                job.many_few_sdatas.many_sdatas.struct.s_vols(i);
            job.subjc(i).maps.mp_vols = {};
            for k = 1:numel(job.many_few_sdatas.many_sdatas.mp_vols)
                job.subjc(i).maps.mp_vols{end+1} = ...
                    job.many_few_sdatas.many_sdatas.mp_vols{k}{i};
            end
        end
    end
end
end
%_______________________________________________________________________

function fout = insert_pref(f,p)
fout = strcat(spm_str_manip(f,'h'),filesep,p,spm_str_manip(f,'t'));
end
%_______________________________________________________________________