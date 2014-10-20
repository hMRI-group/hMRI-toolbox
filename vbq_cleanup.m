function vbq_cleanup(d)
    x=dir(d);
    for i=3:numel(x)
        if x(i).isdir
            if exist(fullfile(d, x(i).name, '_finished_'), 'file')
                disp(['Skipping finished directory ' fullfile(d, x(i).name)]);
            else
                vbq_cleanup(fullfile(d, x(i).name));
            end
        else
            [~,~,ext] = fileparts(x(i).name);
            if strcmp(ext, '.nii')
                idx = strfind(x(i).name, 'in_');
                if ~isempty(idx) && idx(1) == 1 && ~isempty(strfind(x(i).name, '_in.nii'))
                    % Input file, leave alone
                else
                    disp(x(i).name);
                    delete(fullfile(d, x(i).name));
                end
            end
        end
    end
end