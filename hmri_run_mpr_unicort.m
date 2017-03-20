function out = hmri_run_mpr_unicort(job)
% function P = hmri_run_mpr_unicort(P_PDw, P_R1)
% P_PDw: proton density weighted FLASH image (small flip angle image) for
% masking
% P_R1: R1 (=1/T1) map estimated from dual flip angle FLASH experiment
%
% P_R1_unicort: filename of corrected R1 map (same as R1 map with "mh" prefix)
% P_B1: filename of UNICORT estimated B1 map (same as R1 map with "B1_" prefix)
%
% Applies UNICORT correction for RF transmit inhomogoeneities to R1 maps
% estimated from dual angle FLASH experiments. The correction is primarily
% based on the SPM8 "New Segment" toolbox.
% Corrected image and B1+ map is written to the same directory where R1 map is located.
%
% Note: Correction is optimized for 3T Magnetom Tim Trio (Siemens
% Healthcare, Erlangen, Germany) and may need to be re-optimized for other
% field strengths and RF coils. It is was also validated for 1mm isotropic
% whole brain human R1 data only. Smaller coverage and lower resolution may
% lead to suboptimal results. It is recommended to cross-validate with an established B1+
% mapping method (e.g., Lutti et al., MRM 2010) when first applied to different datasets.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
% For theory and validation, see
% Weiskopf et al. (2010), "Unified segmentation based correction of R1
% brain maps for RF transmit field inhomogeneities (UNICORT)", Neuroimage.
%
% Author: N. Weiskopf, WTCN, London
% 29 November 2010

%%

job=hmri_process_data_spec(job);

out.R1 = {};
out.R1u = {};
out.R2s = {};
out.A = {};
out.MT = {};
out.T1w = {};

json = hmri_get_defaults('json');

for ip=1:numel(job.subj)
    P_mtw    = char(job.subj(ip).raw_mpm.MT);
    P_pdw    = char(job.subj(ip).raw_mpm.PD);
    P_t1w    = char(job.subj(ip).raw_mpm.T1);
    
    % determine output directory path
    try 
        cwd = job.subj.output.outdir{1}; % case outdir
    catch 
        Pin = char(job.subj.raw_mpm.MT);
        cwd = fileparts(Pin(1,:)); % case indir
    end
    % save outpath as default for this job
    hmri_get_defaults('outdir',cwd);
    
    [fR1, fR2s, fMT, fA, PPDw, PT1w]  = hmri_MTProt(P_mtw, P_pdw, P_t1w); 
 
    % Retrieve the protocol tag
    hdr = get_metadata(fR1);
    if ~isempty(hdr{1})
        prot_tag = hdr{1}.history.procstep.procpar;
    else
        prot_tag = 'Unknown';
    end
    % Set the tag for the MPMacq set.
    hmri_get_defaults('MPMacq.tag',prot_tag);
    
    % Use default parameters of SPM8 "New Segment" toolbox except for
    % adapted regularization and smoothness of bias field
    % as determined for 3T Magnetom Tim Trio (Siemens Healthcare, Erlangen, Germany)
    % see Weiskopf et al., Neuroimage 2010
    
    
    reg = hmri_get_defaults('unicort.reg');
    FWHM = hmri_get_defaults('unicort.FWHM');
    thr_factor = hmri_get_defaults('unicort.thr');
    
    P_R1     = fR1;
    P_PDw    = PPDw;
    
    % if nargin < 2
    %     P_R1 =[];
    % end
    % if nargin < 1
    %     P_PDw = [];
    % end
    % if isempty(P_PDw)
    %     P_PDw = spm_select(1,'image','Select proton density weighted image');
    % end
    % if isempty(P_R1)
    %     P_R1 = spm_select(1,'image','Select R1 map');
    % end
    
    % create head mask
    V_PDw = spm_vol(P_PDw);
    Y_PDw = spm_read_vols(V_PDw);
    thresh = thr_factor*mode(round(Y_PDw(:))); 
    
    % mask R1 map with head/neck mask
    V_R1 = spm_vol(P_R1);
    Y_R1 = spm_read_vols(V_R1);
    Y_R1 = Y_R1.*(Y_PDw > thresh);
    V_R1_mask = V_R1;
    [p,n,e] = fileparts(V_R1_mask.fname);
    P_R1_mask = fullfile(p,['h' n e]);
    V_R1_mask.fname = P_R1_mask;
    V_R1.descrip = 'Masked R1 map';
    V_R1_mask = spm_write_vol(V_R1_mask,Y_R1);
    
    Vtemp = cat(1,V_PDw,V_R1);
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'map creation';
    Output_hdr.history.procstep.procpar = prot_tag;
    for ctr = 1:numel(Vtemp)
        Output_hdr.history.input{1}.filename = Vtemp(ctr).fname;
        input_hdr = get_metadata(Vtemp(ctr).fname);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{1}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{1}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'Masked R1 map [1000/s]';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(P_R1_mask,Output_hdr,json);
    
    
    %% preparation of spm structure for "New Segment" tool
    
    % clear('matlabbatch');
    tpm_nam = fullfile(spm('dir'),'tpm','enhanced_TPM.nii'); % instead of TPM.nii
    % see http://www.unil.ch/lren/home/menuinst/data--utilities.html
    % Lorio S, Fresard S, Adaszewski S, Kherif F, Chowdhury R, Frackowiak RS, 
    % Ashburner J, Helms G, Weiskopf N, Lutti A, Draganski B. New tissue priors 
    % for improved automated classification of subcortical brain structures on MRI. 
    % Neuroimage. 2016 Apr 15;130:157-66. doi: 10.1016/j.neuroimage.2016.01.062
    
    preproc8.channel.write = [1 1];
    preproc8.tissue(1).tpm = {[tpm_nam ',1']};
    preproc8.tissue(1).ngaus = 2;
    preproc8.tissue(1).native = [0 0];
    preproc8.tissue(1).warped = [0 0];
    preproc8.tissue(2).tpm = {[tpm_nam ',2']};
    preproc8.tissue(2).ngaus = 2;
    preproc8.tissue(2).native = [0 0];
    preproc8.tissue(2).warped = [0 0];
    preproc8.tissue(3).tpm = {[tpm_nam ',3']};
    preproc8.tissue(3).ngaus = 2;
    preproc8.tissue(3).native = [0 0];
    preproc8.tissue(3).warped = [0 0];
    preproc8.tissue(4).tpm = {[tpm_nam ',4']};
    preproc8.tissue(4).ngaus = 3;
    preproc8.tissue(4).native = [0 0];
    preproc8.tissue(4).warped = [0 0];
    preproc8.tissue(5).tpm = {[tpm_nam ',5']};
    preproc8.tissue(5).ngaus = 4;
    preproc8.tissue(5).native = [0 0];
    preproc8.tissue(5).warped = [0 0];
    preproc8.tissue(6).tpm = {[tpm_nam ',6']};
    preproc8.tissue(6).ngaus = 2;
    preproc8.tissue(6).native = [0 0];
    preproc8.tissue(6).warped = [0 0];
    preproc8.warp.reg = [0   0.001   0.5   0.05   0.2];
    preproc8.warp.affreg = 'mni';
    preproc8.warp.samp = 3;
    preproc8.warp.write = [0 0];
    preproc8.warp.mrf = 0;
    % set parameters different from defaults
    preproc8.channel.biasfwhm = FWHM;
    preproc8.channel.biasreg = reg;
    preproc8.channel.vols = {P_R1_mask};
    
    %% run prepared "New Segment" job
    spm_preproc_run(preproc8)
    % spm_jobman('run', matlabbatch);
    clear('matlabbatch');
    
    %% calculate B1+ map from bias field
    [p,n,e] = fileparts(V_R1_mask.fname); %#ok<*NASGU>
    if isempty(spm_select('FPList',p,['^BiasField_' n '.nii']))
        P_biasmap = spm_select('FPList',p,['^BiasField_' n '.img']);
    else
        P_biasmap = spm_select('FPList',p,['^BiasField_' n '.nii']);
    end
    
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'map creation';
    Output_hdr.history.procstep.procpar = prot_tag;
    Output_hdr.history.input{ctr}.filename = P_R1_mask;
    input_hdr = get_metadata(P_R1_mask);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
    Output_hdr.history.output.imtype = 'Bias corrected R1 UNICORT map';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(P_biasmap,Output_hdr,json);
    
    %% create B1+ map from bias field
    V_biasmap = spm_vol(P_biasmap);
    Y_biasmap = spm_read_vols(V_biasmap);
    Y_B1 = sqrt(Y_biasmap)*100.*(Y_PDw > thresh);
    V_B1 = V_R1;
    [p,n,e] = fileparts(V_R1.fname);
    P_B1 = fullfile(p,['B1_' n e]);
    V_B1.fname = P_B1;
    V_B1.descrip = 'B1+ map (p.u. nominal fa)';
    V_B1 = spm_write_vol(V_B1,Y_B1);

    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'map creation';
    Output_hdr.history.procstep.procpar = prot_tag;
    Output_hdr.history.input{1}.filename = V_biasmap.fname;
    input_hdr = get_metadata(V_biasmap.fname);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{1}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{1}.history = 'No history available.';
    end
    Output_hdr.history.output.imtype = 'B1+ map';
    Output_hdr.history.output.units = 'p.u. nominal FA';
    set_metadata(P_B1,Output_hdr,json);
    
    [p,n,e] = fileparts(P_R1_mask);
    P_R1_unicort = fullfile(p, ['m' n e]);
    
    Vtemp = cat(1,V_PDw,V_R1);
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'map creation';
    Output_hdr.history.procstep.procpar = prot_tag;
    for ctr = 1:numel(Vtemp)
        Output_hdr.history.input{ctr}.filename = Vtemp(ctr).fname;
        input_hdr = get_metadata(Vtemp(ctr).fname);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'Bias corrected R1 UNICORT map';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(P_R1_unicort,Output_hdr,json);
    
    out.subj(ip).R1={fullfile(cwd,spm_str_manip(fR1,'t'))};
    out.subj(ip).R1u={fullfile(cwd,spm_str_manip(P_R1_unicort,'t'))};
    out.subj(ip).R2s={fullfile(cwd,spm_str_manip(fR2s,'t'))};
    out.subj(ip).MT={fullfile(cwd,spm_str_manip(fMT,'t'))};
    out.subj(ip).A={fullfile(cwd,spm_str_manip(fA,'t'))};
    out.subj(ip).T1w={fullfile(cwd,spm_str_manip(PT1w,'t'))};
    
    out.R1{end+1} = out.subj(ip).R1{1};
    out.R1u{end+1} = out.subj(ip).R1u{1};
    out.R2s{end+1} = out.subj(ip).R2s{1};
    out.MT{end+1} = out.subj(ip).MT{1};
    out.A{end+1} = out.subj(ip).A{1};
    out.T1w{end+1} = out.subj(ip).T1w{1};
    
    f = fopen(fullfile(cwd, '_finished_'), 'wb');
    fclose(f);
    
end

end