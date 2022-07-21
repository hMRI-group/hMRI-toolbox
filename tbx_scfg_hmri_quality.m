function quality = tbx_scfg_hmri_quality
%==========================================================================
% Configuration (sub)file for the hMRI-toolbox
% PURPOSE
% Implementation of a series of (basic) quality assessment tools,
% mainly based on visual inspection of results, as a first step. 
%
% LIST OF TOOLS
% - check reorientation: display (in CheckReg) the reoriented image
% (preferably the one used for reorientation) and the template used for
% reorientation -> saved as PNG file for quick visual inspection.
% - display all maps: with standard widowing, and saved as PNG file for
% quick wisual inspection.
% - display input images: per image type (e.g. all B1 input images, of all
% FLASH images, etc...), adjust windowing to visualise all images with
% identical windowing. Visual inspection, search for obvious motion
% artefacts, ...
%
% WARNING
% So far, these tools are not meant to automatically flag anomalies. The
% assessment still relies on the expertise of the person visualizing the
% data! Moreover, the PNG files saved don't show everything, and anomalies
% can still keep unnoticed...    
%==========================================================================
% Written by Evelyne Balteau
% Cyclotron Research Centre, University of Liege
%==========================================================================

%--------------------------------------------------------------------------
% res_dir Output directory
%--------------------------------------------------------------------------
res_dir         = cfg_files;
res_dir.tag     = 'res_dir';
res_dir.name    = 'Output directory';
res_dir.val{1}  = {''};
res_dir.help    = {['Files produced by this function will be written into ' ...
    'this output directory. If no directory is given, results images will ' ...
    'be written to the current working directory.']};
res_dir.filter  = 'dir';
res_dir.ufilter = '.*';
res_dir.num     = [0 1];

% ---------------------------------------------------------------------
% check input images
% ---------------------------------------------------------------------
input_images           = cfg_files;
input_images.tag       = 'input_images';
input_images.name      = 'Input images';
input_images.help      = {'Input images.'
    'Select images according to the description of the specific quality check you are about to run. Maximum number of images = 15.'};
input_images.filter    = 'image';
input_images.ufilter   = '.*';
input_images.num       = [0 15];
input_images.val       = {''};

% ---------------------------------------------------------------------
% check_type - Menu listing the types of visual check available 
% ---------------------------------------------------------------------
check_type         = cfg_menu;
check_type.tag     = 'check_type';
check_type.name    = 'Quality check';
check_type.help    = {'Specify the type of quality check you want to run. '
    'You can select either:'
    ['- Input images: Select images all from a single image type ' ...
    '(e.g. all B1 input images, or all spoiled gradient echo (FLASH) images, etc...). ' ...
    'The images will be displayed and the intensity scaling will be ' ...
    'automatically adjusted to visualise all images with ' ...
    'identical intensity scaling. Purpose: visual inspection, search for obvious motion ' ...
    'artefacts, ...']
    ['- Check (re)orientation: Select images that are expected to be realigned ' ...
    '(e.g. FLASH image and template used for auto-reorientation). ' ...
    'Images are displayed with intensity scaling adjusted for each image.']
    ['- qMRI maps: Select B1+ maps, RF sensitivity maps, R1, R2*, MT or PD maps. ' ...
    'The intensity scaling will automatically adjust to use standard and ' ...
    'comparable scaling for each qMRI map. ']};
check_type.labels  = {'Input images'
                'Check (re)orientation'
                'qMRI maps'}';
check_type.values  = {'check_inputs'
                'check_reorient'
                'check_maps'}';
check_type.val     = {'check_inputs'};

% ---------------------------------------------------------------------
% visual_check - Visual inspection of various type of data 
% ---------------------------------------------------------------------
visual_check         = cfg_exbranch;
visual_check.tag     = 'visual_check';
visual_check.name    = 'Visual check';
visual_check.val     = { check_type input_images res_dir };
visual_check.help    = {'Visual inspection of the data.'
    ['Select type of data, images (according to data type) and output directory ' ...
    'where the results will be saved. A PNG image is saved for later review.']};
visual_check.prog    = @hmri_quality_visual_check;
% visual_check.vout    = @vout_quality_visual_check; % not implemented

% ---------------------------------------------------------------------
% quality - list of quality tools
% ---------------------------------------------------------------------
quality         = cfg_choice;
quality.tag     = 'quality';
quality.name    = 'Quality tools';
quality.help    = {
    ['Implementation of a series of (basic) quality assessment tools, ' ...
    'mainly based on visual inspection of results, as a first step. ']};
quality.values  = {visual_check};

end


%==========================================================================
% For all visual check cases...
%==========================================================================
function hmri_quality_visual_check(job)
% job is a structure with fields:
%       check_type: 'check_inputs'/'check_reorient'/'check_maps'
%     input_images: {filenames}
%          res_dir: {dirname}

% flags for hmri_log
log_flags = struct('ComWin',true,'PopUp',false,'LogFile',struct('Enabled',false));

if isempty(job.res_dir)
    job.res_dir = {pwd};
end
job.res_dir = char(job.res_dir);

for cim = 1:length(job.input_images)
    Y{cim} = spm_read_vols(spm_vol(char(job.input_images{cim}))); %#ok<*AGROW>
end

switch(job.check_type)
    case 'check_inputs'
        % define min and max threshold based on percentiles of the ensemble
        % of input images for identical intensity scaling for all images
        tmp = [];
        for cim = 1:length(Y)
            tmp = [tmp; Y{cim}(:)]; 
        end
        mini = prctile(tmp,5);
        maxi = prctile(tmp,95);
%        threshold = (myprctile(mask(:),98)-myprctile(mask(:),2))*0.02+myprctile(mask(:),2);
        for cim=1:length(Y)
            disp_list(cim).fnam = char(job.input_images{cim});
            disp_list(cim).title = spm_file(disp_list(cim).fnam,'basename');
            disp_list(cim).title(strcmp(disp_list(cim).title, '_')) = ' ';
            disp_list(cim).range = [mini maxi];
        end
            
    case 'check_reorient'
        % define min and max threshold based on percentiles for each image
        % separately 
        for cim = 1:length(Y)
            mini = prctile(Y{cim}(:),5);
            maxi = prctile(Y{cim}(:),95);
            disp_list(cim).fnam = char(job.input_images{cim});
            disp_list(cim).title = spm_file(disp_list(cim).fnam,'basename');
            disp_list(cim).title(strcmp(disp_list(cim).title, '_')) = ' ';
            disp_list(cim).range = [mini maxi];
        end
        
    case 'check_maps'
        % define min and max threshold based on data type:
        % - B1: [75-125] p.u.
        % - RF sensitivity: [0 1]
        % - R2*: [0 70] s-1
        % - R1: [0 1.4] s-1
        % - MT saturation: [0 2] p.u.
        % - PD: [50 120] p.u.
        for cim = 1:length(Y)
            notok = false;
            disp_list(cim).fnam = char(job.input_images{cim});
            tmphdr = get_metadata(disp_list(cim).fnam);
            if isempty(tmphdr)
                hmri_log(sprintf(['WARNING (hmri_quality): No metadata associated with qMRI map.' ...
                    '\n%s \nImage displayed with automatic (non-standard) intensity scaling.'], ...
                    disp_list(cim).fnam), log_flags);
                notok = true;
            end
            try
                tmp = tmphdr{1}.history.output.imtype;
                if iscell(tmp)
                    disp_list(cim).title = tmphdr{1}.history.output.imtype{1};
                else
                    disp_list(cim).title = tmphdr{1}.history.output.imtype;
                end                    
            catch %#ok<CTCH>
                hmri_log(sprintf(['WARNING (hmri_quality): No imtype defined for metadata associated with qMRI map.' ...
                    '\n%s \nImage displayed with automatic (non-standard) intensity scaling.'], ...
                    disp_list(cim).fnam), log_flags);
                notok = true;
            end
            if notok 
                mini = prctile(Y{cim}(:),5);
                maxi = prctile(Y{cim}(:),95);
                disp_list(cim).title = spm_file(disp_list(cim).fnam,'basename');
                disp_list(cim).title(strcmp(disp_list(cim).title, '_')) = ' ';
                disp_list(cim).range = [mini maxi];
            else
                if strfind(disp_list(cim).title,'MT') % strfind is case sensitive!
                    disp_list(cim).range = [0 2];
                elseif strfind(disp_list(cim).title,'PD')
                    disp_list(cim).range = [50 120];
                elseif strfind(disp_list(cim).title,'R1')
                    disp_list(cim).range = [0 1.4];
                elseif strfind(disp_list(cim).title,'R2')
                    disp_list(cim).range = [0 70];
                elseif strfind(disp_list(cim).title,'B1')
                    disp_list(cim).range = [75 125];
                elseif strfind(disp_list(cim).title,'RF sensitivity map')
                    disp_list(cim).range = [0 2];
                else % includes case of A maps
                    hmri_log(sprintf(['WARNING (hmri_quality): Non standard scaling for qMRI map.' ...
                        '\n%s \nImage displayed with automatic (non-standard) intensity scaling.'], ...
                        disp_list(cim).fnam), log_flags);
                    mini = prctile(Y{cim}(:),5);
                    maxi = prctile(Y{cim}(:),95);
                    disp_list(cim).range = [mini maxi];
                end
            end         
        end
        
    otherwise
        error('Unknown visual check');
end
        
if ~isempty(disp_list)
    h = hmri_quality_display(disp_list);
    
    colormap(h,'gray');
    spm_orthviews('Redraw');
    res = 1;
    savenam = fullfile(job.res_dir, [job.check_type sprintf('_%0.3d.png', res)]);
    while exist(spm_file(savenam,'suffix','_gray'),'file')
        res = res+1;
        savenam = fullfile(job.res_dir, [job.check_type sprintf('_%0.3d.png', res)]);
    end
    print(h,'-dpng',spm_file(savenam,'suffix','_gray'));
    
    colormap(h,'jet');
    spm_orthviews('Redraw');
    print(h,'-dpng',spm_file(savenam,'suffix','_jet'));
end

end