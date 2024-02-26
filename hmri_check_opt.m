function [opts,ind_repl] = hmri_check_opt(opts_def,opts,loc_opt)

% FORMAT opts = crc_check_flag(opts_def,opts)
%
% Function to automatically check the content of a "opts" structure, using
% a "default opts structure", adding the missing fields and putting in the 
% default value if none was provided.
%
% INPUT:
% opts_def  : default or reference structure
% opts      : input option structure that need to be filled for missing
%             fields with default values
%
% OUPUT:
% opts      : filled up option structure
%_______________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Local option structure:
% loc_opt
%   .verbose    : print out information [true] or not [false, def.] on the 
%                 fixed/updated structure.

if nargin<3
    loc_opt.verbose = false;
end

f_names = fieldnames(opts_def);
% list fields in default structure

Nfields = length(f_names);
ind_repl = zeros(1,Nfields);
for ii=1:Nfields
    % Update the output if
    % - a field is missing
    % - the field is empty when it shouldn't
    if ~isfield(opts,f_names{ii}) || ...
            ( isempty(opts.(f_names{ii})) && ~isempty(opts_def.(f_names{ii})) )
        opts.(f_names{ii}) = opts_def.(f_names{ii});
        ind_repl(ii) = 1;
        if loc_opt.verbose
            fprintf('\n\tAdding field ''%s'' to structure ''%s''.', ...
                f_names{ii},inputname(2));
            jump_line = true;
        end
    end
end

if loc_opt.verbose && jump_line
    fprintf('\n');
end

end