function P_trans = hmri_run_b1map(jobsubj)

%% Processing of B1 maps for B1 bias correction
% FORMAT P_trans = hmri_run_b1map(jobsubj)
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

% retrieve effective acquisition & processing parameters, alternatively
% use defaults 
b1map_params = get_b1map_params(jobsubj);

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
        P = char(jobsubj.raw_fld.b1);
        P_trans  = P(1:2,:);
        
    otherwise 
        fprintf('\nWARNING: unknown B1 type, no B1 map calculation performed.\n');
       
end

% save these into Results directory (nii & json!)
% NOTES: 
%   - if "cleanup" set to true, the B1mapCal directory is deleted when the
%   Map Calculation completes...  
%   - just in case no json files have been saved with the output, the
%   copyfile is called in "try" mode...
%   - must strip the ',1' (at the end of the file extension '.nii,1')
%   otherwise copyfile does not find the files!! 
if ~isempty(P_trans)
    P_trans = spm_file(P_trans,'number','');
    copyfile(P_trans(1,:),fullfile(jobsubj.path.respath, spm_file(P_trans(1,:), 'filename')));
    try copyfile([spm_str_manip(P_trans(1,:),'r') '.json'],fullfile(jobsubj.path.respath, [spm_file(P_trans(1,:), 'basename') '.json'])); end %#ok<*TRYNC>
    copyfile(P_trans(2,:),fullfile(jobsubj.path.respath, spm_file(P_trans(2,:), 'filename')));
    try copyfile([spm_str_manip(P_trans(2,:),'r') '.json'],fullfile(jobsubj.path.respath, [spm_file(P_trans(2,:), 'basename') '.json'])); end
end

% save b1map_params as json-file
spm_jsonwrite(fullfile(jobsubj.path.respath,'MPM_map_creation_b1map_params.json'),b1map_params,struct('indent','\t'));


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
P = char(jobsubj.raw_fld.b1);   % 2 or 4 images
fileTR1 = P(2,:);
fileTR2 = P(1,:);
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

sname = spm_file(V1.fname,'filename');

% save output images
VB1 = V1;
% VB1.pinfo = [max(B1map(:))/16384;0;0];
% VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
% spm_write_vol(VB1,B1map);

% VB1.pinfo = [max(B1map_norm(:))/16384;0;0];
% VB1.fname = fullfile(outpath, [sname '_B1map_norm.nii']);
% spm_write_vol(VB1,B1map_norm);

VB1.pinfo = [max(smB1map_norm(:))/16384;0;0];
VB1.fname = fullfile(outpath, [sname '_smB1map_norm.nii']);
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = cat(1,V2,V1);
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.output.imtype = 'B1+ map (AFI protocol)';
set_metadata(VB1.fname,Output_hdr,json);

% requires anatomic image + map
P_trans  = char(char(fileTR1),char(VB1.fname));

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
% This macro calls the functions hmri_B1Map_unwarp and hmri_B1Map_process
% for correction of image distortions, padding and smoothing of the images.
% Output:
%     - distorted B1 (B1map_*) and error (SDmap_*) maps
%     - undistorted B1 (uB1map_*) and error (uSDmap_*) maps
%     - undistorted, masked and padded B1 maps (muB1map_*)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_*) i.e. FULLY PROCESSED
% At each voxel, this macro selects the 5 pairs of (SE,STE image) (out of
% 11) with maximum signal amplitude in the SE images.
% The sum of square image of all SE images is created (SumOfSq) and
% undistorted (uSumOfSq) for coregistration of the B1 map to an anatomical
% dataset.

json = hmri_get_defaults('json');

P    = char(jobsubj.raw_fld.b1); % B1 data - 11 pairs
Q    = char(jobsubj.raw_fld.b0); % B0 data - 3 volumes

V = spm_vol(P);
n = numel(V);
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
input_files = jobsubj.raw_fld.b1;
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = 'B1+ map calculation (EPI SE/STE protocol)';
 
% save B1 map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'Distorted B1+ map';
Output_hdr.history.output.units = 'p.u.';
V_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','B1 map [%]');
[~,name,e] = fileparts(V_save.fname);
V_save.fname = fullfile(outpath,['B1map_' name e]);
V_save = spm_write_vol(V_save,Y_ab*100);
set_metadata(V_save.fname,Output_hdr,json);

% save SD map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'Distorted SD (error) map';
Output_hdr.history.output.units = 'p.u.';
W_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SD [%]');
W_save.fname = fullfile(outpath,['SDmap_' name e]);
W_save = spm_write_vol(W_save,Y_cd*100);
set_metadata(W_save.fname,Output_hdr,json);

% save SD map (still distorted and not smoothed)
Output_hdr.history.output.imtype = 'SSQ image';
Output_hdr.history.output.units = 'a.u.';
X_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SE SSQ matrix');
X_save.fname = fullfile(outpath,['SumOfSq' name e]);
X_save = spm_write_vol(X_save,Ssq_matrix); %#ok<*NASGU>
set_metadata(X_save.fname,Output_hdr,json);


%-B0 undistortion
%-----------------------------------------------------------------------
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
[fmap_img,unwarp_img] = hmri_B1Map_unwarp(fmfnam, anatfnam, otherfnam, b1map_params);
uanat_img{1} = unwarp_img{1}.fname;
ub1_img{1} = unwarp_img{2}.fname;
ustd_img{1} = unwarp_img{3}.fname;

% set metadata for unwrapped output images
% define generic header for B0-unwarp process
scphasefnam = fullfile(b1map_params.outpath, spm_file(spm_file(fmfnam(2,:),'prefix','sc'),'filename'));
input_files = cat(1,{anatfnam},{fmfnam(1,:)},{fmfnam(2,:)},otherfnam{1},otherfnam{2});
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = 'Unwarp B1 map (EPI SE/STE protocol)';

% set metadata for unwarped B1 image 
Output_hdr.history.output.imtype = 'Unwarped B1 map';
Output_hdr.history.output.units = 'p.u.';
set_metadata(ub1_img{1},Output_hdr,json);

% set metadata for unwarped SD map 
Output_hdr.history.output.imtype = 'Unwarped SD (error) map';
Output_hdr.history.output.units = 'p.u.';
set_metadata(ustd_img{1},Output_hdr,json);

% set metadata for unwarped SSQ map 
Output_hdr.history.output.imtype = 'Unwarped SSQ image';
Output_hdr.history.output.units = 'a.u.';
set_metadata(uanat_img{1},Output_hdr,json);

% set metadata for phase-unwrapped regularised field map (Hz) (fpm_* file) 
Output_hdr.history.output.imtype = 'Phase-unwrapped regularised field map';
Output_hdr.history.output.units = 'Hz';
set_metadata(fmap_img{1}.fname,Output_hdr,json);

% set metadata for Voxel Displacement Map (vdm5_* file) 
Output_hdr.history.output.imtype = 'Voxel displacement map';
Output_hdr.history.output.units = 'Vx';
set_metadata(fmap_img{2}.fname,Output_hdr,json);

% set metadata for phase map scaled between +/-pi (sc* file)
Output_hdr.history.output.imtype = 'Phase map rescaled between [-pi, pi]';
Output_hdr.history.output.units = 'Radians';
set_metadata(scphasefnam,Output_hdr,json);

%-B1 map processing (masking, padding, smoothing, ...)
%--------------------------------------------------------------------------
fpm_img{1} = fmap_img{1};
vdm_img{1} = fmap_img{2};
[allub1_img] = hmri_B1Map_process(ub1_img,ustd_img,vdm_img,fpm_img,b1map_params);

% set metadata for processing B1 images
% define generic header for B1 process
input_files = cat(1,ub1_img,ustd_img,vdm_img{1}.fname,fpm_img{1}.fname);
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = 'Process B1 map (EPI SE/STE protocol)';

% set metadata for each output
for i=1:length(allub1_img)
    Output_hdr.history.output.imtype = allub1_img{i}.descrip;
    Output_hdr.history.output.units = 'p.u.';
    set_metadata(allub1_img{i}.fname,Output_hdr,json);
end

% set correct output for the current subfunction (unwrapped "anatomical"
% image (SSQ) for coregistration and final B1 map):
P_trans  = char(char(uanat_img),char(allub1_img{2}.fname));

end

%% =======================================================================%
% B1 map calculation - SIEMENS tfl_b1map protocol
% Written by Tobias Leutritz (based on calc_AFI_b1map by TL)
%=========================================================================%
function P_trans = calc_tfl_b1map(jobsubj, b1map_params)

json = hmri_get_defaults('json');

VV = char(jobsubj.raw_fld.b1);
P = VV(2,:); % scaled FA map from tfl_b1map sequence
Q = VV(1,:); % FLASH like anatomical from tfl_b1map sequence

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
VB1.fname = fullfile(outpath, [sname '_smB1map_norm.nii']);
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = cat(1,{V2.fname},{V1.fname});
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = 'B1+ map calculation (SIEMENS tfl_b1map protocol)';
set_metadata(VB1.fname,Output_hdr,json);

% requires anatomic image + map
P_trans  = char(char(V2.fname),char(VB1.fname));

end

%% =======================================================================%
% B1 map calculation - SIEMENS rf_map protocol
% Written by Tobias Leutritz 
%=========================================================================%
function P_trans = calc_rf_map(jobsubj, b1map_params)

json = hmri_get_defaults('json');

VV = char(jobsubj.raw_fld.b1);
P = VV(2,:); % scaled FA map from rf_map sequence
Q = VV(1,:); % anatomical image from rf_map sequence

% read header information and volumes
V1 = spm_vol(P); % image volume information
V2 = spm_vol(Q);
Vol1 = spm_read_vols(V1);
Vol2 = spm_read_vols(V2);

% generating the map
B1map_norm = (abs(Vol1)-2048)*180*100/(90*2048); % *100/90 to get p.u.
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
VB1.fname = fullfile(outpath, [sname '_smB1map_norm.nii']);
spm_write_vol(VB1,smB1map_norm);

% set and write metadata
input_files = cat(1,{V2.fname},{V1.fname});
Output_hdr = init_b1_output_metadata(input_files, b1map_params);
Output_hdr.history.procstep.descrip = 'B1+ map calculation (SIEMENS rf_map protocol)';
set_metadata(VB1.fname,Output_hdr,json);

% requires anatomic image + map
P_trans  = char(char(V2.fname),char(VB1.fname));

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
b1_protocol = jobsubj.b1_type;

% define local defaults for b1map
b1map_params.b1type = hmri_get_defaults(['b1map.' b1_protocol '.b1type']); 

% check for existing b1 data
if isempty(jobsubj.raw_fld.b1)
    b1map_params.b1avail = false; % nothing else to be done but send warning if data were expected
    switch(b1map_params.b1type)
        case 'UNICORT'
            b1map_params.procreq = true; % b1 bias correction required
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n'...
                'No B1 map available. UNICORT will be applied.\n']);
        case 'no_B1_correction'
            b1map_params.procreq = false; % no b1 bias correction applied
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                'No B1 map available. No B1 correction applied (semi-quantitative maps only) -----\n']);
        otherwise
            dbstack
            b1map_params.procreq = false;
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                'B1 map calculation cannot proceed because no B1 data available.\n' ...
                'No B1 bias correction will be applied. Results will be semi-\n' ...
                'quantitative only. If you meant to apply B1 bias correction, \n' ...
                'check your data and re-run the batch.']);
    end
else
    switch(b1map_params.b1type)
        case 'UNICORT'
            b1map_params.b1avail = false;
            b1map_params.procreq = true; % b1 bias correction required
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                'B1 map available but UNICORT has been selected and will\n' ...
                'be applied, regardless to the existing B1 data.\n']);
        case 'no_B1_correction'
            b1map_params.b1avail = false;
            b1map_params.procreq = false; % no b1 bias correction applied
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                'B1 map available but no B1 correction applied in agreement \n' ...
                'with the selected B1 processing option. Only semi-quantitative \n' ...
                'maps will be generated.\n']);
        case 'pre_processed_B1'
            b1map_params.b1avail   = true;
            b1map_params.procreq   = false;
            fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                'Preprocessed B1 map available. Assuming it is in percent units. No calculation required.\n']);
        otherwise
            b1map_params.b1avail = true;
            b1map_params.procreq = true;
            
            % retrieve metadata if available
            hdr = get_metadata(jobsubj.raw_fld.b1{1});
            try
                switch(b1map_params.b1type)
                    case 'i3D_EPI'
                        fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                            'SE/STE EPI protocol selected ...\n']);
                        b1map_params.b1acq.beta = get_metadata_val(hdr{1},'B1mapNominalFAValues');
                        b1map_params.b1acq.TM = get_metadata_val(hdr{1},'B1mapMixingTime');
                        b1map_params.b1acq.tert = get_metadata_val(hdr{1},'epiReadoutDuration'); % must take into account PAT but not PF acceleration
                        b1map_params.b1acq.blipDIR = get_metadata_val(hdr{1},'PhaseEncodingDirectionSign');
                        b1map_params.b1proc = hmri_get_defaults(['b1map.' b1_protocol '.b1proc']);
                        % B0 data are required, let's check:
                        if isempty(jobsubj.raw_fld.b0)
                            b1map_params.b0avail = false;
                            fprintf(1,['WARNING: expected B0 map not available for EPI undistortion.\n' ...
                                'No fieldmap correction will be applied.']);
                        else
                            % note that the current implementation assumes that
                            % b0 input images = 2 magnitude images (1st and 2nd
                            % echoes) and 1 presubtracted phase image.
                            b1map_params.b0avail = true;
                            b1map_params.b0acq.shortTE = get_metadata_val(jobsubj.raw_fld.b0{1},'EchoTime');
                            b1map_params.b0acq.longTE = get_metadata_val(jobsubj.raw_fld.b0{2},'EchoTime');
                            b1map_params.b0acq.iformat = 'PM';
                        end
                        
                    case 'i3D_AFI'
                        fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                            'AFI protocol selected ...\n']);
                        tr = get_metadata_val(hdr{1},'RepetitionTimes');
                        b1map_params.b1acq.TR2TR1ratio = tr(2)/tr(1);
                        b1map_params.b1acq.alphanom = get_metadata_val(hdr{1},'FlipAngle');
                        
                    case 'tfl_b1_map'
                        fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                            'SIEMENS tfl_b1map protocol selected ...\n']);
                        
                    case 'rf_map'
                        fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                            'SIEMENS rf_map protocol selected ...\n']);
                end
            catch
                fprintf(1, ['\n---------------- B1 MAP CALCULATION ----------------\n' ...
                    'WARNING: possibly no metadata associated to the input images. \n' ...
                    'Default acquisition and processing parameters will be used.\n' ...
                    '%s data type is assumed with protocol %s.\n'], b1map_params.b1type, b1_protocol);
                b1map_def = hmri_get_defaults(['b1map.' b1_protocol]);
                f = fieldnames(b1map_def);
                for cfi=1:length(f)
                    if ~isfield(b1map_params,f{cfi})
                        b1map_params.(f{cfi}) = b1map_def.(f{cfi});
                    end
                end
            end
            fprintf(1,'\n');
            if isfield(b1map_params, 'b1acq')
                fprintf(1,'B1 acquisition parameters:\n');
                disp(b1map_params.b1acq);
            end
            if isfield(b1map_params, 'b0acq')
                fprintf(1,'B0 acquisition parameters:\n');
                disp(b1map_params.b0acq);
            end
            if isfield(b1map_params, 'b1proc')
                fprintf(1,'B1 processing parameters:\n');
                disp(b1map_params.b1proc);
            end
    end
        
end
end

%=========================================================================%
% To arrange the metadata structure for B1 map calculation output.
%=========================================================================%
function metastruc = init_b1_output_metadata(input_files, b1map_params)

proc.descrip = 'B1+ map calculation';
proc.version = hmri_get_version;
proc.params = b1map_params;

output.imtype = 'B1+ map';
output.units = 'p.u.';

metastruc = init_output_metadata_structure(input_files, proc, output);

end