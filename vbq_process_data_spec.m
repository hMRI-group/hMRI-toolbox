function job = vbq_process_data_spec(job)
if isfield(job.data_spec, 'sdata_multi')
    % ---- MT ----
    s = 1;
    job.subj(s).output = job.data_spec.sdata_multi.output;
    job.subj(s).raw_mpm.MT = {};
    for i=1:numel(job.data_spec.sdata_multi.raw_mpm.MT)
        fname = job.data_spec.sdata_multi.raw_mpm.MT{i};
        [~, name, ext] = fileparts(fname);
        
        if isempty(fname) || strcmp([name ext], 'separator.nii')
            s = s + 1;
            job.subj(s).output = job.data_spec.sdata_multi.output;
            job.subj(s).raw_mpm.MT = {};
        else
            job.subj(s).raw_mpm.MT{end+1} = fname;
        end
    end
    % ---- PD ----
    s = 1;
    job.subj(s).raw_mpm.PD = {};
    for i=1:numel(job.data_spec.sdata_multi.raw_mpm.PD)
        fname = job.data_spec.sdata_multi.raw_mpm.PD{i};
        [~, name, ext] = fileparts(fname);
        
        if isempty(fname) || strcmp([name ext], 'separator.nii')
            s = s + 1;
            job.subj(s).raw_mpm.PD = {};
        else
            job.subj(s).raw_mpm.PD{end+1} = fname;
        end
    end
    % ---- T1 ----
    s = 1;
    job.subj(s).raw_mpm.T1 = {};
    for i=1:numel(job.data_spec.sdata_multi.raw_mpm.T1)
        fname = job.data_spec.sdata_multi.raw_mpm.T1{i};
        [~, name, ext] = fileparts(fname);
        
        if isempty(fname) || strcmp([name ext], 'separator.nii')
            s = s + 1;
            job.subj(s).raw_mpm.T1 = {};
        else
            job.subj(s).raw_mpm.T1{end+1} = fname;
        end
    end
    % ---- raw_fld ----
    if isfield(job.data_spec.sdata_multi, 'raw_fld')
        % ---- b0 ----
        s = 1;
        job.subj(s).raw_fld.b0 = {};
        for i=1:numel(job.data_spec.sdata_multi.raw_fld.b0)
            fname = job.data_spec.sdata_multi.raw_fld.b0{i};
            [~, name, ext] = fileparts(fname);

            if isempty(fname) || strcmp([name ext], 'separator.nii')
                s = s + 1;
                job.subj(s).raw_fld.b0 = {};
            else
                job.subj(s).raw_fld.b0{end+1} = fname;
            end
        end
        % ---- b1 ----
        s = 1;
        job.subj(s).raw_fld.b1 = {};
        for i=1:numel(job.data_spec.sdata_multi.raw_fld.b1)
            fname = job.data_spec.sdata_multi.raw_fld.b1{i};
            [~, name, ext] = fileparts(fname);

            if isempty(fname) || strcmp([name ext], 'separator.nii')
                s = s + 1;
                job.subj(s).raw_fld.b1 = {};
            else
                job.subj(s).raw_fld.b1{end+1} = fname;
            end
        end
    end
    
    if isempty(job.subj(end).raw_mpm.MT)
        job.subj = job.subj(1:(end-1));
    end
else
    job.subj = job.data_spec.subj;
end
