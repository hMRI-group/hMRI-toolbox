function P_trans = hmri_create_b1map(jobsubj)

%% Processing of B1 maps for B1 bias correction
% FORMAT P_trans = hmri_create_b1map(jobsubj)
%    jobsubj - are parameters for one subject out of the job list.
%    NB: ONE SINGLE DATA SET FROM ONE SINGLE SUBJECT IS PROCESSED HERE,
%    LOOP OVER SUBJECTS DONE AT HIGHER LEVEL.
%    P_trans - a vector of file names with P_trans(1,:) = anatomical volume
%        for coregistration and P_trans(2,:) = B1 map in percent units.
%_______________________________________________________________________
% Written by E. Balteau, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
%_______________________________________________________________________
% Modified by T. Leutritz in 2016 in order to use the SIEMENS product
% sequences 'rf_map' and 'tfl_b1map'. The latter produces essentially
% a FLASH like image and a flip angle map (multiplied by 10) based on
% Chung S. et al.: "Rapid B1+ Mapping Using a Preconditioning RF Pulse with
% TurboFLASH Readout", MRM 64:439-446 (2010).
%_______________________________________________________________________

flags = jobsubj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ CREATE B1 MAP - %s.m (%s) ============', mfilename, datestr(now)),flags);

% retrieve effective acquisition & processing parameters, alternatively
% use defaults 
b1map_params = get_b1map_params(jobsubj);

% save b1map_params as json-file
spm_jsonwrite(fullfile(jobsubj.path.supplpath,'hMRI_map_creation_b1map_params.json'),b1map_params,struct('indent','\t'));

% init output
P_trans = [];

% return if nothing else to be done (no B1 correction or UNICORT cases)
if ~b1map_params.b1avail
    return;
end

% calculate B1 map according to b1 data type
switch(b1map_params.b1type)
    case 'i3D_AFI'
        % processing B1 map from AFI data
        P_trans  = calc_AFI_b1map(jobsubj, b1map_params);
        
    case 'i3D_EPI'
        % processing B1 map from SE/STE EPI data
        P_trans  = calc_SESTE_b1map(jobsubj, b1map_params);
        
    case 'tfl_b1_map'
        % processing B1 map from tfl_b1map data
        P_trans  = calc_tfl_b1map(jobsubj, b1map_params);
        
    case 'rf_map'
        % processing B1 map from rf_map data
        P_trans  = calc_rf_map(jobsubj, b1map_params);
        
    case 'pre_processed_B1'
        if b1map_params.scafac ~= 1
            % rescale if scaling factor other than 1 is provided
            rescaled_fnam = fullfile(jobsubj.path.b1path, spm_file(spm_file(b1map_params.b1input(2,:),'suffix','_rescaled'),'filename'));
            calcflags.descrip = sprintf('Pre-processed B1 map rescaled with factor %f', b1map_params.scafac);
            outcalc = spm_imcalc(b1map_params.b1input(2,:),rescaled_fnam,sprintf('%f*i1',b1map_params.scafac),calcflags);
            % set and write metadata
            json = hmri_get_defaults('json');
            input_files = b1map_params.b1input(2,:);
            Output_hdr = init_b1_output_metadata(input_files, b1map_params);
            Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (Rescaling)'];
            Output_hdr.history.output.imtype = sprintf('Pre-processed B1 map rescaled with factor %f', b1map_params.scafac);
            set_metadata(rescaled_fnam,Output_hdr,json);
            % replace original B1 map by the rescaled one
            b1map_params.b1input = char(b1map_params.b1input(1,:), rescaled_fnam);
        end
        P_trans  = b1map_params.b1input(1:2,:);
        
    otherwise 
        hmri_log(sprintf('WARNING: unknown B1 type, no B1 map calculation performed.'),b1map_params.defflags);
       
end

% copy P_trans output to Results/Supplementary directory (nii & json!) and
% make P_trans point to the copied files (so coregistration is applied to
% them).
%
% NOTES: 
%   - if "cleanup" set to true, the B1mapCalc directory is deleted when the
%   Map Calculation completes...  
%   - just in case no json files have been saved with the output, the
%   copyfile is called in "try" mode...
%   - must strip the ',1' (at the end of the file extension '.nii,1')
%   otherwise copyfile does not find the files!! 

if ~isempty(P_trans)
    P_trans = spm_file(P_trans,'number','');
    P_trans_copy{1} = fullfile(jobsubj.path.b1respath, spm_file(P_trans(1,:), 'filename'));
    P_trans_copy{2} = fullfile(jobsubj.path.b1respath, spm_file(P_trans(2,:), 'filename'));
    copyfile(deblank(P_trans(1,:)), P_trans_copy{1});
    try copyfile([spm_str_manip(P_trans(1,:),'r') '.json'],[spm_str_manip(P_trans_copy{1},'r') '.json']); end %#ok<*TRYNC>
    copyfile(deblank(P_trans(2,:)), P_trans_copy{2});
    try copyfile([spm_str_manip(P_trans(2,:),'r') '.json'],[spm_str_manip(P_trans_copy{2},'r') '.json']); end
    P_trans = char(P_trans_copy{1},P_trans_copy{2});
end

hmri_log(sprintf('\t============ CREATE B1 MAP: completed (%s) ============', datestr(now)),b1map_params.nopuflags);

end

%% =======================================================================%
% B1 map calculation - AFI protocol
%=========================================================================%
function P_trans = calc_AFI_b1map(jobsubj, b1map_params)

% default format specifications for the output metadata
json = hmri_get_defaults('json');

% define output dir
outpath = jobsubj.path.b1path;
b1map_params.outpath = outpath;

% NB: both phase and magnitude images can be provided but only the
% magnitude images (first series) are used. Phase images (second series)
% are not used. In each series, first image = TR2 (long TR) and second
% image = TR1 (short TR).
fileTR1 = b1map_params.b1input(2,:);
fileTR2 = b1map_params.b1input(1,:);
V1 = spm_vol(fileTR1);
V2 = spm_vol(fileTR2);
Y1 = spm_read_vols(V1);
Y2 = spm_read_vols(V2);

TR1 = 1; % only the ratio [TR2/TR1=n] matters
TR2 = b1map_params.b1acq.TR2TR1ratio;
alphanom = b1map_params.b1acq.alphanom;

% Mask = squeeze(Vol1);
% threshold = (prctile(Mask(:),98)-prctile(Mask(:),2))*0.1+prctile(Mask(:),2);
% Mask = (Mask>threshold);

B1map = acos((Y2./Y1*TR2/TR1-1)./(TR2/TR1*ones(size(Y1))-Y2./Y1))*180/pi;
B1map_norm = abs(B1map)*100/alphanom;

% smoothed map
smB1map_norm = zeros(size(B1map_norm));
pxs = sqrt(sum(V1.mat(1:3,1:3).^2)); % Voxel resolution
smth = 8./pxs;
spm_smooth(B1map_norm,smB1map_norm,smth);

% masking
% B1map = B1map.*Mask;
% B1map_norm = B1map_norm.*Mask;
% smB1map_norm = smB1map_norm.*Mask;

sname = spm_file(V1.fname,'basename');

% save output images
VB1 = V1;
% VB1.pinfo = [max(B1map(:))/16384;0;0];
% VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
% spm_write_vol(VB1,B1map);

VB1.pinfo = [max(B1map_norm(:))/16384;0;0];
VB1.descrip = 'B1+ map - normalised (p.u.) - AFI protocol';
VB1.fname = fullfile(outpath, [sname '_unsmoothed_B1map.nii']);
spm_write_vol(VB1,B1map_norm);

VB1.pinfo = [max(smB1map_norm(:))/16384;0;0];
VB1.descrip = 'B1+ map - smoothed and normalised (p.u.) - AFI protocol';
VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = b1map_params.b1input;
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (AFI protocol)'];
Output_hdr.history.output.imtype = 'B1+ map (AFI protocol)';
set_metadata(VB1.fname,Output_hdr,json);

% Rename anatomical reference for uniformity between protocols
B1ref = fullfile(outpath, [sname '_B1ref.nii']);
copyfile(char(fileTR1),B1ref);
try copyfile([spm_str_manip(char(fileTR1),'r') '.json'],[spm_str_manip(B1ref,'r') '.json']); end %#ok<*TRYNC>

% requires anatomic image + map
P_trans  = char(B1ref,char(VB1.fname));

% VB1.fname = fullfile(outpath, [sname '_B1map_mask.nii']);
% spm_write_vol(VB1,Mask);

end

%% =======================================================================%
% B1 map calculation - SE/STE EPI protocol
%=========================================================================%
function P_trans = calc_SESTE_b1map(jobsubj, b1map_params)
% Calculation of B1 maps based on 3D EPI spin echo (SE) and stimulated
% (STE) echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping
% Input: 11 pairs of (SE, STE) images for B1 map calculation and 3 images
% for B0 map calculation.
% This macro calls the functions hmri_create_B1Map_unwarp and
% hmri_create_B1Map_process for correction of image distortions, padding
% and smoothing of the images. 
% Output:
%     - distorted B1 (B1map_*) and error (SDmap_*) maps
%     - undistorted B1 (uB1map_*) and error (uSDmap_*) maps
%     - undistorted, masked and padded B1 maps (muB1map_*)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_*)
%                                                   i.e. FULLY PROCESSED
% At each voxel, this macro selects the 5 pairs of (SE,STE image) (out of
% 11) with maximum signal amplitude in the SE images.
% The sum of square image of all SE images is created (SumOfSq) and
% undistorted (uSumOfSq) for coregistration of the B1 map to an anatomical
% dataset.
% 
% For coherence among B1 protocols, the fully processed B1 map (smuB1map_*)
% is renamed *_B1map.nii, while the undistorted SoS image (uSumOfSq) is
% renamed *_B1ref for anatomical reference.

json = hmri_get_defaults('json');

P    = b1map_params.b1input; % B1 data - 11 pairs
Q    = b1map_params.b0input; % B0 data - 3 volumes

V = spm_vol(P);
n = numel(V);

assert(rem(n,2)==0, ...
   ['B1 mapping image volumes must be a set of SE, STE pairs ' ...
   'thus the number of input volumes (currently %d) must be even.'], n);

assert(n == 2 * numel(b1map_params.b1acq.beta), ...
   ['Number of B1 mapping image pairs (%d) does not match ' ...
    'the number of nominal flip angles (%d)!'], ...
    n/2, numel(b1map_params.b1acq.beta));

% splitting images into SE and STE volumes
% TODO: check that correct data selected (compare echo times?)
V_SE = V(1:2:end);
V_STE = V(2:2:end);

% calc_SESTE_b1map expects fa in decreasing order
[b1map_params.b1acq.beta, fa_order] = sort(b1map_params.b1acq.beta, 'descend');

% rearranging volumes in decreasing fa
V_SE = V_SE(fa_order);
V_STE = V_STE(fa_order);

% TODO: Keeping SE and STE separated might be easier to read
% but requires a lot of code changes below, so merging into
% original list with corrected order
V(1:2:end) = V_SE;
V(2:2:end) = V_STE;

Y_tmptmp = zeros([V(1).dim(1:2) n]);
Y_ab = zeros(V(1).dim(1:3));
Y_cd = zeros(V(1).dim(1:3));
Index_Matrix = zeros([V(1).dim(1:3) b1map_params.b1proc.Nonominalvalues]);
real_Y_tmp = zeros([V(1).dim(1:2) 2*b1map_params.b1proc.Nonominalvalues]);

Ssq_matrix = sqrt(sum(spm_read_vols(V(1:2:end)).^2,4));

%-Define output directory
%-----------------------------------------------------------------------
outpath = jobsubj.path.b1path;
b1map_params.outpath = outpath;

%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',V(1).dim(3),'B1 map fit','planes completed');

%-Loop over planes computing result Y
%-----------------------------------------------------------------------
clear Temp_mat;
corr_fact = exp(b1map_params.b1acq.TM/b1map_params.b1proc.T1);
for p = 1:V(1).dim(3) %loop over the partition dimension of the data set
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
    for i = 1:n/2
        M = inv(B*inv(V(1).mat)*V(1).mat); %#ok<*MINV>
        Y_tmptmp(:,:,((i-1)*2+1))  = real( ...
            acos(corr_fact*spm_slice_vol(V((i-1)*2+2),M,V(1).dim(1:2),0) ./ ...
            (spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0)+b1map_params.b1proc.eps))/pi*180/b1map_params.b1acq.beta(i) ...
            ); % nearest neighbor interpolation
        Y_tmptmp(:,:,((i-1)*2+2))  = 180/b1map_params.b1acq.beta(i) - Y_tmptmp(:,:,((i-1)*2+1));
        Temp_mat(:,:,i) = spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0); %#ok<*AGROW>
    end
    
    [~,indexes] = sort(Temp_mat,3);
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
            for k=1:b1map_params.b1proc.Nonominalvalues
                real_Y_tmp(x_nr,y_nr,2*k-1) = Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1)-1);
                real_Y_tmp(x_nr,y_nr,2*k)   = Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1));
                Index_Matrix(x_nr,y_nr,p,k) = indexes(x_nr,y_nr,indexes(x_nr,y_nr,n/2-k+1));
            end
        end
    end
    
    Y_tmp = sort(real(real_Y_tmp), 3); % take the real value due to noise problems
    Y_sd  = zeros([V(1).dim(1:2) (b1map_params.b1proc.Nonominalvalues+1)]);
    Y_mn  = zeros([V(1).dim(1:2) (b1map_params.b1proc.Nonominalvalues+1)]);
    for i = 1:(b1map_params.b1proc.Nonominalvalues+1)
        Y_sd(:,:,i) = std(Y_tmp(:,:,i:(i + b1map_params.b1proc.Nonominalvalues-1)), [], 3);
        Y_mn(:,:,i) = mean(Y_tmp(:,:,i:(i + b1map_params.b1proc.Nonominalvalues-1)), 3);
    end
    
    [~,min_index] = min(Y_sd,[],3); % !! min_index is a 2D array. Size given by resolution along read and phase directions
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
            Y_ab(x_nr,y_nr,p) = Y_mn(x_nr,y_nr, min_index(x_nr,y_nr)); % Y_ab is the relative flip angle value averaged over the n flip angles (determined by minizing the SD i.e. keeping the most uniform relative flip angle values)
            Y_cd(x_nr,y_nr,p) = Y_sd(x_nr,y_nr, min_index(x_nr,y_nr)); % Y_cd is the corresponding standard deviation between the relative flip angle values
        end
    end
    spm_progress_bar('Set',p);
end

%-Save everything in OUTPUT dir
%-----------------------------------------------------------------------
% define generic output header
input_files = b1map_params.b1input;
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (EPI SE/STE protocol)'];
 
% save B1 map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Distorted B1+ map';
Output_hdr.history.output.units = 'p.u.';
V_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','B1 map [%]');
[~,outname,e] = fileparts(V_save.fname);
V_save.fname = fullfile(outpath,['B1map_' outname e]);
V_save = spm_write_vol(V_save,Y_ab*100);
set_metadata(V_save.fname,Output_hdr,json);

% save SD map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Distorted SD (error) map';
Output_hdr.history.output.units = 'p.u.';
W_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SD [%]');
W_save.fname = fullfile(outpath,['SDmap_' outname e]);
W_save = spm_write_vol(W_save,Y_cd*100);
set_metadata(W_save.fname,Output_hdr,json);

% save SD map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - SSQ image';
Output_hdr.history.output.units = 'a.u.';
X_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SE SSQ matrix');
X_save.fname = fullfile(outpath,['SumOfSq' outname e]);
X_save = spm_write_vol(X_save,Ssq_matrix); %#ok<*NASGU>
set_metadata(X_save.fname,Output_hdr,json);


%-B0 undistortion
%-----------------------------------------------------------------------
% since B0 data will be coregistered and resliced with the B1 data, we copy
% them into the calcpath directory to avoid altering the the raw data:
Qtmp = cell(size(Q,1),1);
for i=1:size(Q,1)
    Qtmp{i} = fullfile(outpath, spm_file(Q(i,:), 'filename'));
    copyfile(deblank(Q(i,:)), Qtmp{i});
    try copyfile([spm_str_manip(deblank(Q(i,:)),'r') '.json'],[spm_str_manip(Qtmp{i},'r') '.json']); end %#ok<*TRYNC>
end
Q = char(Qtmp);

% magnitude image 
% NOTE: must strip the ',1' (at the end of the file extension '.nii,1')!!
magfnam = spm_file(Q(1,:),'number','');
% phase image
phasefnam = spm_file(Q(3,:),'number','');
% both fieldmap images
fmfnam = char(phasefnam,magfnam);
% image to be corrected ("anatomical" reference = SSQ image)
anatfnam = X_save.fname;
% other images to be corrected (distorted B1 and SD maps)
otherfnam{1} = V_save.fname;
otherfnam{2} = W_save.fname;

% unwarp
[fmap_img,unwarp_img] = hmri_create_B1Map_unwarp(fmfnam, anatfnam, otherfnam, b1map_params);
uanat_img{1} = unwarp_img{1}.fname;
ub1_img{1} = unwarp_img{2}.fname;
ustd_img{1} = unwarp_img{3}.fname;

% set metadata for unwrapped output images
% define generic header for B0-unwarp process
scphasefnam = fullfile(b1map_params.outpath, spm_file(spm_file(fmfnam(2,:),'prefix','sc'),'filename'));
% relate outputs to inputs remaining visible after cleanup! i.e. original
% B1 and B0 mapping images (not to the intermediate images created
% during B1 calculation):
% input_files = cat(1,{anatfnam},{fmfnam(1,:)},{fmfnam(2,:)},otherfnam{1},otherfnam{2});
input_files = char(b1map_params.b1input,b1map_params.b0input);
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (EPI SE/STE protocol)'];

% set metadata for unwarped B1 image 
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Unwarped B1 map';
Output_hdr.history.output.units = 'p.u.';
set_metadata(ub1_img{1},Output_hdr,json);

% set metadata for unwarped SD map 
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Unwarped SD (error) map';
Output_hdr.history.output.units = 'p.u.';
set_metadata(ustd_img{1},Output_hdr,json);

% set metadata for unwarped SSQ map 
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Unwarped SSQ image for anatomical reference';
Output_hdr.history.output.units = 'a.u.';
set_metadata(uanat_img{1},Output_hdr,json);

% set metadata for phase-unwrapped regularised field map (Hz) (fpm_* file) 
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Phase-unwrapped regularised field map';
Output_hdr.history.output.units = 'Hz';
set_metadata(fmap_img{1}.fname,Output_hdr,json);

% set metadata for Voxel Displacement Map (vdm5_* file) 
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Voxel displacement map';
Output_hdr.history.output.units = 'Vx';
set_metadata(fmap_img{2}.fname,Output_hdr,json);

% set metadata for phase map scaled between +/-pi (sc* file)
Output_hdr.history.output.imtype = 'SE/STE B1 mapping - Phase map rescaled between [-pi, pi]';
Output_hdr.history.output.units = 'Radians';
set_metadata(scphasefnam,Output_hdr,json);

%-B1 map processing (masking, padding, smoothing, ...)
%--------------------------------------------------------------------------
fpm_img{1} = fmap_img{1};
vdm_img{1} = fmap_img{2};
[allub1_img] = hmri_create_B1Map_process(ub1_img,ustd_img,vdm_img,fpm_img,b1map_params);

% set metadata for processing B1 images
% define generic header for B1 process
% relate outputs to inputs remaining visible after cleanup! i.e. original
% B1 and B0 mapping images (not to the intermediate images created
% during B1 calculation):
% input_files = cat(1,ub1_img,ustd_img,vdm_img{1}.fname,fpm_img{1}.fname);
input_files = char(b1map_params.b1input,b1map_params.b0input);
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (EPI SE/STE protocol)'];

% set metadata for each output
for i=1:length(allub1_img)
    Output_hdr.history.output.imtype = ['SE/STE B1 mapping - ' allub1_img{i}.descrip];
    Output_hdr.history.output.units = 'p.u.';
    set_metadata(allub1_img{i}.fname,Output_hdr,json);
end

% set correct output for the current subfunction (unwrapped "anatomical"
% image (SSQ) for coregistration and final B1 map). For coherence among B1
% protocol, rename these files *_B1ref (for anatomical reference) and
% *_B1map (for B1+ bias map in p.u.):  
B1map = fullfile(outpath,[outname '_B1map.nii']);
copyfile(allub1_img{2}.fname, B1map);
try copyfile([spm_str_manip(allub1_img{2}.fname,'r') '.json'],[spm_str_manip(B1map,'r') '.json']); end
B1ref = fullfile(outpath,[outname '_B1ref.nii']);
copyfile(uanat_img{1}, B1ref);
try copyfile([spm_str_manip(uanat_img{1},'r') '.json'],[spm_str_manip(B1ref,'r') '.json']); end
P_trans  = char(B1ref, B1map);

end

%% =======================================================================%
% B1 map calculation - SIEMENS tfl_b1map protocol
% Written by Tobias Leutritz (based on calc_AFI_b1map by TL)
%=========================================================================%
function P_trans = calc_tfl_b1map(jobsubj, b1map_params)

json = hmri_get_defaults('json');

P = b1map_params.b1input(2,:); % scaled FA map from tfl_b1map sequence
Q = b1map_params.b1input(1,:); % anatomical image from tfl_b1map sequence

% read header information and volumes
V1 = spm_vol(P); % image volume information
V2 = spm_vol(Q);
Vol1 = spm_read_vols(V1);
Vol2 = spm_read_vols(V2);

alphanom = get_metadata_val(P,'FlipAngle'); % nominal flip angle of tfl_b1map

% generating the map
B1map_norm = abs(Vol1)*10/alphanom;

% smoothed map
smB1map_norm = zeros(size(B1map_norm));
pxs = sqrt(sum(V1.mat(1:3,1:3).^2)); % Voxel resolution
smth = 8./pxs;
spm_smooth(B1map_norm,smB1map_norm,smth);

% Save everything in OUTPUT dir
%-----------------------------------------------------------------------
% determine output directory path
outpath = jobsubj.path.b1path;
b1map_params.outpath = outpath;

sname = spm_file(V1.fname,'basename');

VB1 = V1;
VB1.pinfo = [max(smB1map_norm(:))/16384;0;0]; % what is this for? (TL)
VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
VB1.descrip = 'Smoothed & normalised (p.u.) B1 bias map - TFL B1map protocol';
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = cat(1,{V2.fname},{V1.fname});
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (SIEMENS tfl_b1map protocol)'];

set_metadata(VB1.fname,Output_hdr,json);

% copy also anatomical image to outpath to prevent modification of original data
anat_fname = fullfile(outpath, [spm_file(V2.fname, 'basename') '_B1ref.nii']);
copyfile(V2.fname, anat_fname);
try copyfile([spm_str_manip(V2.fname,'r') '.json'],[spm_str_manip(anat_fname,'r') '.json']); end %#ok<*TRYNC>

% requires anatomic image + map
P_trans  = char(char(anat_fname),char(VB1.fname));

end

%% =======================================================================%
% B1 map calculation - SIEMENS rf_map protocol
% Written by Tobias Leutritz 
%=========================================================================%
function P_trans = calc_rf_map(jobsubj, b1map_params)

json = hmri_get_defaults('json');

P = b1map_params.b1input(2,:); % scaled FA map from rf_map sequence
Q = b1map_params.b1input(1,:); % anatomical image from rf_map sequence

% read header information and volumes
V1 = spm_vol(P); % image volume information
V2 = spm_vol(Q);
Vol1 = spm_read_vols(V1);
Vol2 = spm_read_vols(V2);
alphanom = get_metadata_val(P,'FlipAngle'); % nominal flip angle of rf_map

% generating the map
B1map_norm = (abs(Vol1)-2048)*180*100/(alphanom*2048); % *100/alpha to get p.u.
% the formula (abs(Vol1)-2048)*180/2048 would result in an absolute FA map

% smoothed map
smB1map_norm = zeros(size(B1map_norm));
pxs = sqrt(sum(V1.mat(1:3,1:3).^2)); % Voxel resolution
smth = 8./pxs;
spm_smooth(B1map_norm,smB1map_norm,smth);

% Save everything in OUTPUT dir
%-----------------------------------------------------------------------
% determine output directory path
outpath = jobsubj.path.b1path;
b1map_params.outpath = outpath;

sname = spm_file(V1.fname,'basename');

VB1 = V1;
VB1.pinfo = [max(smB1map_norm(:))/16384;0;0]; % what is this for? (TL)
VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
VB1.descrip = 'Smoothed & normalised (p.u.) B1 bias map - TFL B1map protocol';
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = cat(1,{V2.fname},{V1.fname});
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (SIEMENS rf_map protocol)'];
set_metadata(VB1.fname,Output_hdr,json);

% copy also anatomical image to outpath to prevent modification of original data
anat_fname = fullfile(outpath, [spm_file(V2.fname, 'basename') '_B1ref.nii']);
copyfile(V2.fname, anat_fname);
try copyfile([spm_str_manip(V2.fname,'r') '.json'],[spm_str_manip(anat_fname,'r') '.json']); end %#ok<*TRYNC>

% requires anatomic image + map
P_trans  = char(char(anat_fname),char(VB1.fname));

end


%% =======================================================================%
% Determine whether b1 data are available and whether any processing should
% be applied. If so, all the required parameters for b1map calculation are
% retrieved, including b1map and b0map acquisition parameters and
% processing parameters, if applicable. Check whether input data are
% coherent with the processing type selected. Missing parameters will be 
% retrieved from the hmri_get_defaults.
%=========================================================================%
function b1map_params = get_b1map_params(jobsubj)

% retrieve b1 protocol from job 
% (can be different - a variation of - the b1 type)
f = fieldnames(jobsubj.b1_type);
b1_protocol = f{1};

% pre-set filename of defaults file
deffnam = '';
custom_def = false;

% load customized defaults parameters from customized defaults file if any
% (the customized defaults file must be run to overwrite the standard
% defaults parameters)
if isfield(jobsubj.b1_type.(b1_protocol),'b1parameters')
    % first reinitialise processing parameters to standard defaults:
    hmri_b1_standard_defaults;
    deffnam = fullfile(fileparts(mfilename('fullpath')),'config','hmri_b1_standard_defaults.m');
    custom_def = false;

    % then, if customized defaults file available, run it to overwrite
    % standard defaults parameters.
    if isfield(jobsubj.b1_type.(b1_protocol).b1parameters,'b1defaults')
        deffnam = jobsubj.b1_type.(b1_protocol).b1parameters.b1defaults;
        spm('Run',deffnam);
        custom_def = true;
    end
end

% load all B1 bias correction defaults parameters & add default file
b1map_params = hmri_get_defaults(['b1map.' b1_protocol]); 
b1map_params.defaults_file = deffnam;
b1map_params.custom_defaults = custom_def;

% flags for logging information and warnings
b1map_params.defflags = jobsubj.log.flags; % default flags
b1map_params.nopuflags = jobsubj.log.flags; % force no Pop-Up
b1map_params.nopuflags.PopUp = false; 

hmri_log(sprintf('\t------------ B1 MAP CALCULATION (%s) %s ------------',b1_protocol, datestr(now)),b1map_params.nopuflags);

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
b1map_params.SPMver = sprintf('%s (%s)', v, r);

% load B1 input images if any
% (NB: if a 'b1input' field is present, it should NOT be empty)
if isfield(jobsubj.b1_type.(b1_protocol),'b1input')
    b1map_params.b1input = char(spm_file(jobsubj.b1_type.(b1_protocol).b1input,'number',''));        
    if isempty(b1map_params.b1input)
        hmri_log(sprintf(['WARNING: expected B1 input images missing. Switching to "no \n' ...
            '\tB1 correction" mode. If you meant to apply B1 bias correction, \n' ...
            '\tcheck your data and re-run the batch.']),b1map_params.defflags);
        b1_protocol = 'no_B1_correction';
        b1map_params = hmri_get_defaults('b1map.no_B1_correction'); 
    end
end
        
% load B0 input images if any
% (NB: if a 'b0input' field is present, it may be empty)
if isfield(jobsubj.b1_type.(b1_protocol),'b0input')
    b1map_params.b0input = char(spm_file(jobsubj.b1_type.(b1_protocol).b0input,'number',''));
    if isempty(b1map_params.b0input)
        % hmri_log(sprintf(['WARNING: expected B0 fieldmap not available for EPI undistortion.\n' ...
        %     '\tNo fieldmap correction will be applied.']),b1map_params.defflags);
        % b1map_params.b0avail = false; 
        hmri_log(sprintf(['WARNING: expected B0 fieldmap not available for EPI undistortion.\n' ...
            '\tThe current implementation does not allow you to apply EPI-based B1 bias \n' ...
            '\tcorrection without phase unwrapping. Switching to "no B1 correction" mode.\n' ...
            '\tIf you meant to apply B1 bias correction, check your data and re-run the batch.']),b1map_params.defflags);
        b1_protocol = 'no_B1_correction';
        b1map_params = hmri_get_defaults('b1map.no_B1_correction');        
    end
end

% process job inputs according to B1 type
switch b1_protocol
    case 'UNICORT'
        hmri_log(sprintf('No B1 map available. UNICORT will be applied.'),b1map_params.nopuflags);

    case 'no_B1_correction'
        hmri_log(sprintf('No B1 map available. No B1 correction applied (semi-quantitative maps only)'),b1map_params.nopuflags);
        
    case 'pre_processed_B1'
        b1map_params.scafac = jobsubj.b1_type.(b1_protocol).scafac;
        if ~isempty(b1map_params.b1input)
            if b1map_params.scafac == 1
                hmri_log(sprintf('Preprocessed B1 map available. \nAssuming it is in percent units of the nominal flip angle. \nNo calculation required.'),b1map_params.defflags);
            else
                hmri_log(sprintf('Preprocessed B1 map available. \nScaling factor provided: %f. Assuming B1 map will be expressed \nin p.u. of the nominal flip angle after rescaling.', b1map_params.scafac),b1map_params.defflags);
            end
        end

    case 'i3D_EPI'
        if ~isempty(b1map_params.b1input)
            hmri_log(sprintf('SE/STE EPI protocol selected ...'),b1map_params.nopuflags);
            b1hdr = get_metadata(b1map_params.b1input(1,:));
            
            try
                tmp = get_metadata_val(b1hdr{1},'B1mapNominalFAValues');
                if isempty(tmp)
                    hmri_log(sprintf('WARNING: using defaults value for nominal SE/STE flip angle values \n(%s) instead of metadata', ...
                        sprintf('%d ',b1map_params.b1acq.beta)),b1map_params.defflags);
                else b1map_params.b1acq.beta = tmp; 
                end
                
                tmp = get_metadata_val(b1hdr{1},'B1mapMixingTime');
                if isempty(tmp)
                    hmri_log(sprintf('WARNING: using defaults value for mixing time \n(%d ms) instead of metadata', ...
                    b1map_params.b1acq.TM),b1map_params.defflags);
                else b1map_params.b1acq.TM = tmp; 
                end
                
                tmp = get_metadata_val(b1hdr{1},'epiReadoutDuration'); % must take into account PAT but not PF acceleration
                if isempty(tmp)
                    hmri_log(sprintf('WARNING: using defaults value for EPI readout duration\n(%d ms) instead of metadata', ...
                        b1map_params.b1acq.tert),b1map_params.defflags);
                else b1map_params.b1acq.tert = tmp; 
                end
                
                tmp = get_metadata_val(b1hdr{1},'PhaseEncodingDirectionSign');
                if isempty(tmp)
                    hmri_log(sprintf('WARNING: using defaults value for PE direction\n(%d) instead of metadata', ...
                        b1map_params.b1acq.blipDIR),b1map_params.defflags);
                else b1map_params.b1acq.blipDIR = tmp; 
                end
                
                % consistency check for T1 value and field strength:
                tmp = get_metadata_val(b1hdr{1},'MagneticFieldStrength');
                supportedB0 = false;
                matchT1fieldstrength = false;
                if ~isempty(tmp)
                    switch round(tmp)
                        case 3
                            supportedB0 = true;
                            expectedT1 = 1192;
                        case 7
                            supportedB0 = true;
                            expectedT1 = 1633;
                        otherwise
                            supportedB0 = false;
                            expectedT1 = NaN;
                    end
                    if b1map_params.b1proc.T1 == expectedT1
                        matchT1fieldstrength = true;
                    end
                    if ~supportedB0
                        hmri_log(sprintf(['WARNING: field strength (B0 = %.0fT) not supported. The reference T1' ...
                            '\nvalue for B1 map calculation for that field strength is not currently ' ...
                            '\nimplemented in the hMRI-toolbox. Please make sure the assumed ' ...
                            '\nvalue (T1 = %.0f ms) is correct, otherwise set it via a customised ' ...
                            '\nB1 default file (config/local/hmri_b1_local_defaults.m).' ...
                            '\nIf the value is already properly set, just ignore this message.'], ...
                            tmp, b1map_params.b1proc.T1),b1map_params.defflags);
                    else
                        if ~matchT1fieldstrength && custom_def
                            hmri_log(sprintf(['WARNING: the assumed T1 value for B1 map calculation does not ' ...
                                '\nmatch the expected value for the used field strength: ' ...
                                '\n    B0 = %.0fT, T1 = %d/%d (expected/actual) ms.' ...
                                '\n\nPlease check T1 value is properly set in your local settings ' ...
                                '\n(see hmri_def.b1map.i3D_EPI.b1proc.T1 in your customised ' ...
                                '\n%s config file).' ...
                                '\n\nRecommended values are: ' ...
                                '\n    - @3T: T1 = 1192 ms' ...
                                '\n    - @7T: T1 = 1633 ms' ...
                                '\n\nIf the value was set differently on purpose, just ignore this message.'], ...
                                tmp, expectedT1, b1map_params.b1proc.T1, char(spm_file(deffnam,'filename'))), b1map_params.defflags);
                        elseif ~matchT1fieldstrength && ~custom_def
                             hmri_log(sprintf(['WARNING: the assumed T1 value for B1 map calculation ' ...
                                '\nhas been set to match the used field strength: ' ...
                                '\n    B0 = %.0fT, T1 = %d ms.' ...
                                '\n\nPlease consider to properly set the T1 value uing a local ' ...
                                '\ndefaults file (see config/local/hmri_b1_local_defaults.m ' ...
                                '\nand parameter hmri_def.b1map.i3D_EPI.b1proc.T1 therein).' ...
                                '\n\nRecommended values are: ' ...
                                '\n    - @3T: T1 = 1192 ms' ...
                                '\n    - @7T: T1 = 1633 ms'], ...
                                tmp, expectedT1),b1map_params.defflags);
                            b1map_params.b1proc.T1 = expectedT1;
                        end
                    end
                    b1map_params.b1proc.matchT1fieldstrength = matchT1fieldstrength;
                    b1map_params.b1proc.expectedT1 = expectedT1;
                end
                    
                
                if ~isempty(b1map_params.b0input)
                    % note that the current implementation assumes that
                    % b0 input images = 2 magnitude images (1st and 2nd
                    % echoes) and 1 presubtracted phase image.
                    tmp = get_metadata_val(b1map_params.b0input(1,:),'EchoTime');
                    if isempty(tmp)
                        hmri_log(sprintf('WARNING: using defaults value for B0 mapping TEs\n(short TE=%.2fms) instead of metadata', ...
                            b1map_params.b0acq.shortTE),b1map_params.defflags);
                    else b1map_params.b0acq.shortTE = tmp; 
                    end
                    
                    tmp = get_metadata_val(b1map_params.b0input(2,:),'EchoTime');
                    if isempty(tmp)
                        hmri_log(sprintf('WARNING: using defaults value for B0 mapping TEs\n(long TE=%.2fms) instead of metadata', ...
                            b1map_params.b0acq.longTE),b1map_params.defflags);
                    else b1map_params.b0acq.longTE = tmp; 
                    end
                    b1map_params.b0acq.iformat = 'PM';
                end
            catch %#ok<*CTCH>
                hmri_log(sprintf(['WARNING: possibly no metadata associated to the input images. \n' ...
                    'Default acquisition and processing parameters will be used.']),b1map_params.defflags);
            end
        end
    case 'i3D_AFI'
        if ~isempty(b1map_params.b1input)
            hmri_log(sprintf('AFI protocol selected ...'),b1map_params.nopuflags);
            b1hdr = get_metadata(b1map_params.b1input(1,:));
            
            try
                tr = get_metadata_val(b1hdr{1},'RepetitionTimes');
                if isempty(tr)
                    hmri_log(sprintf('WARNING: using defaults values for TRs\n(TR ratio = %.1f) instead of metadata', ...
                        b1map_params.b1acq.TR2TR1ratio),b1map_params.defflags);
                else b1map_params.b1acq.TR2TR1ratio = tr(2)/tr(1); 
                end
                
                tmp = get_metadata_val(b1hdr{1},'FlipAngle');
                if isempty(tmp)
                    hmri_log(sprintf('WARNING: using defaults value for flip ange \n(%d deg) instead of metadata', ...
                        b1map_params.b1acq.alphanom), b1map_params.defflags);
                else b1map_params.b1acq.alphanom = tmp; 
                end
            catch
                hmri_log(sprintf(['WARNING: possibly no metadata associated to the input images. \n' ...
                    'Default acquisition and processing parameters will be used.']),b1map_params.defflags);
            end
        end
        
    case 'tfl_b1_map'
        if ~isempty(b1map_params.b1input)
            hmri_log(sprintf('SIEMENS tfl_b1map protocol selected ...'),b1map_params.nopuflags);
        end
                        
    case 'rf_map'
        if ~isempty(b1map_params.b1input)
            hmri_log(sprintf('SIEMENS rf_map protocol selected ...'),b1map_params.nopuflags);
        end
        
    otherwise
        hmri_log(sprintf(['WARNING: something must have gone wrong in the JOB configuration.\n' ...
            '\tUnknown B1 processing methods, assuming "no B1 correction" mode.']),b1map_params.defflags);
        b1_protocol = 'no_B1_correction';
        b1map_params = hmri_get_defaults('b1map.no_B1_correction');
end

% print acquisition and processing parameters
if isfield(b1map_params, 'b1acq')
    hmri_log(sprintf('B1 acquisition parameters (check carefully!):\n\n%s', ...
        printstruct(b1map_params.b1acq)),b1map_params.defflags);
end
if isfield(b1map_params, 'b0acq')
    hmri_log(sprintf('B0 acquisition parameters (check carefully!):\n\n%s', ...
        printstruct(b1map_params.b0acq)),b1map_params.defflags);
end
if isfield(b1map_params, 'b1proc')
    hmri_log(sprintf('B1 processing parameters (check carefully!):\n\n%s', ...
        printstruct(b1map_params.b1proc)),b1map_params.defflags);
end

end

%=========================================================================%
% To arrange the metadata structure for B1 map calculation output.
%=========================================================================%
function metastruc = init_b1_output_metadata(input_files, b1map_params)

proc.descrip = ['hMRI toolbox - ' mfilename '.m - B1+ map calculation'];
proc.version = hmri_get_version;
proc.params = b1map_params;

output.imtype = 'B1+ map';
output.units = 'p.u.';

metastruc = init_output_metadata_structure(input_files, proc, output);

end

%=========================================================================%
% To rpint a structure into text - assumes simple structure (no
% sub-structure in it at this point) 
%=========================================================================%
function s = printstruct(struc)

s = '';
fntmp = fieldnames(struc);
for cf = 1:length(fntmp)
    s = sprintf('%s %16s: %s\n', s, fntmp{cf}, num2str(struc.(fntmp{cf})));
end
end
