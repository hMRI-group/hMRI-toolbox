function out = hmri_create_unicort(P_PDw, P_R1, jobsubj)
% function P = hmri_create_unicort(P_PDw, P_R1, jobsubj)
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
flags.def = jobsubj.log.flags; % default flags
flags.nopu = jobsubj.log.flags; % force no Pop-Up
flags.nopu.PopUp = false; 

hmri_log(sprintf('\t============ UNICORT: CORRECT R1 MAPS FOR B1+ BIAS - %s.m (%s) ============', mfilename, datetime('now')),flags.nopu);

% json metadata default options
json = hmri_get_defaults('json');

% Use default parameters of SPM8 "New Segment" toolbox except for
% adapted regularization and smoothness of bias field
% as determined for 3T Magnetom Tim Trio (Siemens Healthcare, Erlangen, Germany)
% see Weiskopf et al., Neuroimage 2010

% first reinitialise processing parameters to standard defaults:
hmri_b1_standard_defaults;

% if customized defaults file available, run it to overwrite standard
% defaults parameters:
if isfield(jobsubj.b1_type.UNICORT.b1parameters,'b1defaults')
    deffnam = jobsubj.b1_type.UNICORT.b1parameters.b1defaults;
    spm('Run',deffnam);
end
unicort_params = hmri_get_defaults('b1map.UNICORT.procpar');
reg = unicort_params.reg;
FWHM = unicort_params.FWHM;
thr_factor = unicort_params.thr;

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
unicort_params.SPMver = sprintf('%s (%s)', v, r);

% output directories
mpmpath = jobsubj.path.mpmpath;
b1path = jobsubj.path.b1path;
respath = jobsubj.path.respath;
supplpath = jobsubj.path.supplpath;

% create head mask
V_PDw = spm_vol(P_PDw);
Y_PDw = spm_read_vols(V_PDw);
thresh = thr_factor*mode(round(Y_PDw(:)));

% mask R1 map with head/neck mask
V_R1 = spm_vol(P_R1);
Y_R1 = spm_read_vols(V_R1);
Y_R1 = Y_R1.*(Y_PDw > thresh);
V_R1_mask = V_R1;
outfnam = spm_file(P_R1,'filename');
P_R1_mask = fullfile(mpmpath,spm_file(outfnam,'suffix','_masked'));
V_R1_mask.fname = P_R1_mask;
V_R1_mask.descrip = 'Masked R1 map';
spm_write_vol(V_R1_mask,Y_R1);

% set and save metadata
input_files = char(P_PDw,P_R1);
Output_hdr = init_unicort_output_metadata(input_files, unicort_params);
Output_hdr.history.output.imtype = 'Masked R1 map';
Output_hdr.history.output.units = 's-1';
set_metadata(P_R1_mask,Output_hdr,json);


%% preparation of job for Unified Segmentation
% uniform defaults used across the toolbox, in particular, enhanced TPMs
% (provided with the toolbox) are preferred instead of SPM's TPM.nii -
% see http://www.unil.ch/lren/home/menuinst/data--utilities.html
% Lorio S, Fresard S, Adaszewski S, Kherif F, Chowdhury R, Frackowiak RS,
% Ashburner J, Helms G, Weiskopf N, Lutti A, Draganski B. New tissue priors
% for improved automated classification of subcortical brain structures on MRI.
% Neuroimage. 2016 Apr 15;130:157-66. doi: 10.1016/j.neuroimage.2016.01.062

job_US = hmri_get_defaults('segment');
job_US.channel.vols = {P_R1_mask};
job_US.channel.biasreg = reg;
job_US.channel.biasfwhm = FWHM;
job_US.channel.write = [1 1]; % write both BiasField and BiasCorrected volume
for ctis=1:length(job_US.tissue)
    job_US.tissue(ctis).native = [0 0]; % no need to write c* volumes
end
% The following were the parameter used in 2010 for the UNICROT paper, now
% updated to current recommendations for unified segmentation (SPM12 - John
% Ashburner)
% job_US.tissue(1).ngaus = 2; % default is 1
% job_US.tissue(2).ngaus = 2; % default is 1
% job_US.warp.mrf = 0; % default is 1
% job_US.warp.cleanup = 0; % default is 1

%% run prepared "New Segment" job
output_list = spm_preproc_run(job_US);

%% calculate B1+ map from bias field
P_biasmap = output_list.channel.biasfield{1};

% set and save metadata
input_files = P_R1_mask;
Output_hdr = init_unicort_output_metadata(input_files, unicort_params);
Output_hdr.history.output.imtype = 'BiasField for corrected R1 UNICORT map calculation';
Output_hdr.history.output.units = 'a.u.';
set_metadata(P_biasmap,Output_hdr,json);

%% create B1+ map from bias field
V_biasmap = spm_vol(P_biasmap);
Y_biasmap = spm_read_vols(V_biasmap);
Y_B1 = sqrt(Y_biasmap)*100.*(Y_PDw > thresh);
V_B1 = V_R1;
P_B1 = fullfile(b1path,spm_file(outfnam,'suffix','_B1map'));
V_B1.fname = P_B1;
V_B1.descrip = 'UNICORT-estimated B1+ map (p.u.)';
spm_write_vol(V_B1,Y_B1);

% set and save metadata
input_files = P_biasmap;
Output_hdr = init_unicort_output_metadata(input_files, unicort_params);
Output_hdr.history.output.imtype = 'UNICORT-estimated B1+ map';
Output_hdr.history.output.units = 'p.u.';
set_metadata(P_B1,Output_hdr,json);

% Bias corrected R1 map
P_R1_unicort = output_list.channel.biascorr{1};

% set and save metadata
input_files = char(P_PDw,P_R1);
Output_hdr = init_unicort_output_metadata(input_files, unicort_params);
Output_hdr.history.output.imtype = 'R1 map corrected for B1+ bias (UNICORT)';
Output_hdr.history.output.units = 's-1';
set_metadata(P_R1_unicort,Output_hdr,json);

% define output file names
out.R1u = {fullfile(respath,spm_file(outfnam,'suffix','_UNICORT'))};
out.B1u = {fullfile(supplpath,spm_str_manip(P_B1,'t'))};

% now copy files from calc directory into results directory (nii & json!)
% NB: only final maps that will be further analysed are kept in Results.
% Other maps (e.g. uncorrected R1) are moved to Results/Supplemntary. 
% 1. copy UNICORT R1 map into Results 
copyfile(P_R1_unicort,out.R1u{1});
try copyfile([spm_str_manip(P_R1_unicort,'r') '.json'],[spm_str_manip(out.R1u{1},'r') '.json']); end %#ok<*TRYNC>
% 2. move R1 (not bias corrected) into Results/Supplementary
movefile(P_R1,fullfile(supplpath,outfnam));
try movefile([spm_str_manip(P_R1,'r') '.json'],fullfile(supplpath,[spm_file(outfnam,'basename') '.json'])); end %#ok<*TRYNC>
% 3. copy UNICORT-estimated B1 map into Results/Supplementary
copyfile(P_B1,out.B1u{1});
try copyfile([spm_str_manip(P_B1,'r') '.json'],[spm_str_manip(out.B1u{1},'r') '.json']); end %#ok<*TRYNC>

% save unicort params as json-file
spm_jsonwrite(fullfile(supplpath, 'hMRI_map_creation_unicort_params.json'),unicort_params,struct('indent','\t'));

hmri_log(sprintf('\t============ UNICORT: completed (%s) ============', datetime('now')),flags.nopu);

end

%% =======================================================================%
% To arrange the metadata structure for B1 map calculation output.
%=========================================================================%
function metastruc = init_unicort_output_metadata(input_files, unicort_params)

proc.descrip = ['hMRI toolbox - ' mfilename '.m - UNICORT for B1+ bias estimation'];
proc.version = hmri_get_version;
proc.params = unicort_params;

output.imtype = 'B1+ map';
output.units = 'p.u.';

metastruc = init_output_metadata_structure(input_files, proc, output);

end