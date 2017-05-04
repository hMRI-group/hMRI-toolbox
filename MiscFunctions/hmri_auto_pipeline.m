function job=hmri_auto_pipeline(job)
if ~isfield(job.auto_pipeline, 'auto_pipeline_yes')
    return
end

sInDir = char(job.auto_pipeline.auto_pipeline_yes.auto_pipeline_dir);
if sInDir(end) == '\'
    sInDir = sInDir(1:end-1);
end
sOutDir = char(job.auto_pipeline.auto_pipeline_yes.auto_pipeline_odir);
if sOutDir(end) == '\'
    sOutDir = sOutDir(1:end-1);
end
bUnpack = job.auto_pipeline.auto_pipeline_yes.auto_pipeline_unpack;
bCreateHierarchy = job.auto_pipeline.auto_pipeline_yes.auto_pipeline_hierarchy;

if exist(fullfile(sOutDir, 'progress.mat'), 'file')
    progress = load(fullfile(sOutDir, 'progress.mat'));
else
    progress = struct('stage', 1);
end

if bUnpack && progress.stage <= 1
    files = list_files_rec(sInDir);
    for i=1:numel(files)
        if ~isempty(regexp(files{i}, '.tar$', 'match'))
            untar(files{i}, fileparts(files{i}));
            delete(files{i});
        end
    end
end

if bCreateHierarchy && progress.stage <= 2
    %         if progress.stage < 2
    %             progress.stage = 2;
    %             files = list_files_rec(sInDir);
    %             f = fopen(fullfile(sOutDir, 'hmri_files.txt'), 'wt');
    %             for i = 1:numel(files)
    %                 fwrite(f, sprintf('%s\n', files{i}));
    %             end
    %             fclose(f);
    %             save(fullfile(sOutDir, 'progress.mat'), '-struct', 'progress');
    %         end
    
    progress.stage = 2;
    save(fullfile(sOutDir, 'progress.mat'), '-struct', 'progress');
    
    cmd1 = ['java -jar "' fullfile(spm('dir'), 'toolbox', 'hmri', 'Dicomymizer.jar') '" anonymizer -hier PatientName:StudyDate:ProtocolName:SeriesNumber_SeriesDescription -outdir "' sOutDir '" -indir "' sInDir '" -nc -sv'];
    disp(cmd1);
    system(cmd1);
end

progress.stage = 3; %#ok<STRNU>
save(fullfile(sOutDir, 'progress.mat'), '-struct', 'progress');
hmri_cleanup(sOutDir); % remove results from previous run

subj_count = 0;
subj_orig = job.subj(1);
pat = dir(sOutDir);
for i=3:numel(pat)
    if ~pat(i).isdir
        continue
    end
    subj = subj_orig;
    stud = dir(fullfile(sOutDir, pat(i).name));
    for k=3:numel(stud)
        P = fullfile(sOutDir, pat(i).name, stud(k).name);
        seq = dir(P);
        
        for m=3:numel(seq)
            P2 = fullfile(P, seq(m).name);
            if ~should_convert(job, seq(m).name)
                disp(['Skipping ' P2 ' ...']);
                continue
            end
            ser = dir(P2);
            ser = num_sort_dir(ser);
            if ~isempty(regexp(seq(m).name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_b0, 'match'))
                N = numel(ser);
            else
                N = 1;
            end
            for n=1:N % Convert just first series, the one we will use for hMRI
                P3 = fullfile(P2, ser(n).name);
                
                mosaic_result = '';
                if job.auto_pipeline.auto_pipeline_yes.auto_pipeline_mosaic
                    mosaic_result = process_mosaic(sOutDir, P3);
                end
                
                if isempty(mosaic_result)
                    old_dir = pwd;
                    cd(P3);
                    status = local_dicom_convert(P3);
                    if status ~= 0
                        error(['Problem converting DICOM to Nifti : ' num2str(status)]);
                    end
                    all_files = dir();
                    for jj=3:numel(all_files)
                        [~, name, ext] = fileparts(all_files(jj).name);
                        if strcmp(ext, '.nii')
                            if local_is_input_file(all_files(jj).name)
                                continue
                            end
                            fix_origin(all_files(jj).name);
                            movefile(all_files(jj).name, ['in_' name '_in' ext]);
                        else
                            delete(all_files(jj).name);
                        end
                    end
                    cd(old_dir);
                end
            end
        end
        
        seq_names = {seq.name};
        
        subj.raw_mpm.MT = list_files_rec(multi_fullfile(P, find_str(seq_names, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_mt)), 1);
        subj.raw_mpm.PD = list_files_rec(multi_fullfile(P, find_str(seq_names, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_pd)), 1);
        subj.raw_mpm.T1 = list_files_rec(multi_fullfile(P, find_str(seq_names, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_t1)), 1);
        
        finished = false;
        for ii=1:numel(subj.raw_mpm.MT)
            if strcmp(subj.raw_mpm.MT{ii}, '_finished_')
                finished = true;
                break;
            end
        end
        
        if finished
            disp(['Skipping ' pat(i).name ' : maps have already been computed.']);
            break;
        end
        
        check_count('MT', subj.raw_mpm.MT, [6 8]);
        check_count('PD', subj.raw_mpm.PD, 8);
        check_count('T1', subj.raw_mpm.T1, [6 8]);
        
        if isfield(subj, 'raw_fld')
            subj.raw_fld.b1 = list_files_rec(multi_fullfile(P, find_str(seq_names, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_b1)));
            subj.raw_fld.b0 = list_files_rec(multi_fullfile(P, find_str(seq_names, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_b0)));
            
            check_count('B1', subj.raw_fld.b1, 22);
            check_count('B0', subj.raw_fld.b0, 3);
        end
    end
    
    if ~isempty(subj.raw_mpm.MT) && ~isempty(subj.raw_mpm.PD) && ~isempty(subj.raw_mpm.T1) && (~isfield(subj, 'raw_fld') || (~isempty(subj.raw_fld.b1) && ~isempty(subj.raw_fld.b0)))
        subj_count = subj_count + 1;
        job.subj(subj_count) = subj;
    else
        disp(['Missing images for ' pat(i).name]);
    end
end
end

function check_count(name, list, expected)
if numel(expected) == 1
    expected(2) = expected(1);
end
n = numel(list);
if n < expected(1) || n > expected(2)
    error([num2str(n) ' instead of expected ' num2str(expected(1)) ' - ' num2str(expected(2)) ' in ' name]);
end
end

function x = num_sort_dir(x)
x=x(3:end);
ary = zeros(numel(x),1);
for i=1:numel(x)
    tmp = sscanf(x(i).name, '%f', 1);
    if isempty(tmp)
        ary(i) = 1;
    else
        ary(i) = tmp;
    end
end
[~, idx] = sort(ary);
x = x(idx);
end

function F = list_files_rec(D, num_dir)
if ~exist('num_dir', 'var')
    num_dir = Inf;
end
cnt_dir = num_dir;
F = {};
if iscell(D) && numel(D)>1
    for j=1:numel(D)
        F = [F; list_files_rec(D{j}, num_dir)]; %#ok<AGROW>
    end
    return
end
D = char(D);
x = dir(D);
x=x(3:end);
ary = {x.name};
for j=1:numel(ary)
    ary{j} = sscanf(ary{j}, '%f', 1);
    if isempty(ary{j})
        ary{j} = 1;
    end
end
[~, idx] = sort(cell2mat(ary));
x=x(idx);
count = 0;
for j=1:numel(x)
    if ~x(j).isdir
        count = count + 1;
        F{count,1} = fullfile(D, x(j).name); %#ok<AGROW>
    end
end
for j=1:numel(x)
    if x(j).isdir && cnt_dir>0
        cnt_dir = cnt_dir - 1;
        F = [F; list_files_rec(fullfile(D, x(j).name), num_dir)]; %#ok<AGROW>
    end
end
end

function remove_empty_dir(D)
D = char(D);
x = dir(D);
for j=3:numel(x)
    if x(j).isdir
        remove_empty_dir(fullfile(D, x(j).name));
    end
end
try
    rmdir(D)
catch e %#ok<NASGU>
    % leave alone non-empty directories
end
end

function S = find_str(A, R)
R = char(R);
S = {};
count = 0;
for j=1:numel(A)
    if ~isempty(regexp(A{j}, R, 'match'))
        count = count + 1;
        S{count} = A{j}; %#ok<AGROW>
    end
end
S=sort(S);
if ~isempty(S)
    S=S{1};
else
    S='null';
end
end

function M = multi_fullfile(P, F)
M = {};
if ~iscell(F)
    F = {F};
end
for j=1:numel(F)
    M{j} = fullfile(P, F{j}); %#ok<AGROW>
end
end

function res = process_mosaic(sOutDir, path)
res = GetImgFromMosaic(sOutDir, path);
x=dir(fullfile(path, 'Echo*'));
oldwd = pwd;
cd(path);
for j=1:size(x,1)
    p = fullfile(path, x(j).name);
    status1 = local_dicom_convert(p);
    if status1 ~= 0
        error('Problem converting DICOM to Nifti');
    end
    y = dir(p);
    for k=3:numel(y)
        [~,~,ext] = fileparts(y(k).name);
        if ~strcmp(ext, '.nii') % Delete Dicoms
            delete(fullfile(p, y(k).name));
        end
    end
    remove_empty_dir(p);
end
all_nii = dir('*.nii');
for j=1:numel(all_nii)
    if local_is_input_file(all_nii(j).name)
        continue
    end
    fix_origin(all_nii(j).name);
    [~, name, ext] = fileparts(all_nii(j).name);
    movefile(all_nii(j).name, ['in_' name '_in' ext]);
end
cd(oldwd);
end

function status = local_dicom_convert(dir_name)
old_dir = pwd();
cd(dir_name);
x = dir();
files = {};
for i=3:numel(x)
    [~,~,ext] = fileparts(x(i).name);
    if ~strcmp(ext, '.nii') && ~x(i).isdir
        files{end+1} = x(i).name; %#ok<AGROW>
    end
end
hdr = spm_dicom_headers(char(files), true);
cd(old_dir);
out = spm_dicom_convert(hdr, 'all', 'flat', 'nii');
if numel(out.files) > 0
    status = 0;
else
    status = 1;
end
end

function status = local_dicom_convert2(dir_name)
cpu = computer;
if strcmp(cpu, 'PCWIN') || strcmp(cpu, 'PCWIN64')
    exename = 'dcm2nii';
elseif strcmp(cpu, 'GLNX86')
    exename = 'dcm2nii.glnx86';
elseif strcmp(cpu, 'GLNXA64')
    exename = 'dcm2nii.glnxa64';
else
    error('dcm2nii : unsupported architecture');
end
x = dir(dir_name);
all_nii = true;
for i=3:numel(x)
    [~,~,ext] = fileparts(x(i).name);
    if ~strcmp(ext, '.nii')
        all_nii = false;
        break;
    end
end
if all_nii
    disp(['Only Nifti files found. Not doing conversion for ' dir_name]);
    status = 0;
    return
end
exename = fullfile(spm('dir'), 'toolbox', 'hmri', 'dcm2nii', exename);
cmd = [exename ' -b "' fullfile(spm('dir'), 'toolbox', 'hmri', 'dcm2nii', 'dcm2nii.ini') '" -o . -g N "' dir_name '"'];
status = system(cmd);
end

function fix_origin(file)
image2set_hdr=spm_vol(file);
orig_mat = spm_get_space(file);
real_pos1 = orig_mat * [image2set_hdr.dim/2 1]';
mat = spm_matrix([-real_pos1(1), -real_pos1(2), -real_pos1(3), 0, 0, 0, 1, 1, 1]);
spm_get_space(file, mat * orig_mat);
end

function ret = should_convert(job, protocol_name)
ret = false;
if ~isempty(regexp(protocol_name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_mt, 'match')) || ...
        ~isempty(regexp(protocol_name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_pd, 'match')) || ...
        ~isempty(regexp(protocol_name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_t1, 'match')) || ...
        (isfield(job.subj(1), 'raw_fld') && ...
        (~isempty(regexp(protocol_name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_b1, 'match')) || ...
        ~isempty(regexp(protocol_name, job.auto_pipeline.auto_pipeline_yes.auto_pipeline_b0, 'match'))))
    ret = true;
end
end

function ret = local_is_input_file(fname)
ret = false;
[~,name,ext] = fileparts(fname);
if strcmp(ext, '.nii')
    idx = strfind(name, 'in_');
    if ~isempty(idx) && idx(1) == 1 && ~isempty(strfind(fname, '_in.nii'))
        ret = true;
    end
end
end
