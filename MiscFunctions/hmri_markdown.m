function hmri_markdown
% Extract information from the toolbox m-files and output them as text that
% can be copied-pasted into the github wiki (markdown format).
%
% There are 2 types of m2markdown operations:
% 1. converting the job configuration tree, i.e. *_cfg_* and *_scfg_* files
%    defining the batching interface into a series entries and
%    corresponding help sections. 
% 2. converting the help header of the functions into markdown-formatted
%    text. 
%
% The output text files are then copied-pasted into the Help section in the
% wiki. 
%
% File derived from that of the SPM8 distribution.
% http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Turning the cfg files into a text file
if ~nargin,
    % this is just to check whether cfg tools have been already loaded in
    % the working space or need to be loaded before proceeding...
    try
        isa('cfg','cfg_choice'); % this will return 0 if cfg_choice is known, and throw an error otherwise
    catch %#ok<CTCH>
        spm_jobman('initcfg'); % load all cfg tools (not doing it systematically since it is time consuming)
    end
    c = tbx_cfg_hmri;
end

hMRIdir = fileparts(spm_str_manip(mfilename('fullpath'),'h'));

for ii=1:numel(c.values),
    bn = c.values{ii}.tag;
    fp = fopen(fullfile(hMRIdir,'manual',['batch_',bn,'.txt']),'w');
    if fp==-1; return; end;
    chapter(c.values{ii},fp);
    fclose(fp);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. picking all the functions help files and put them into functions.tex
fp = fopen(fullfile(hMRIdir,'manual','functions_help_section.txt'),'w');
if fp==-1, return; end;
excluded = {'MiscFunctions','manual','tbx_cfg','tbx_scfg','etpm'};

% List all scripts
f = spm_select('FPListRec',hMRIdir,'.*\.m$');
f = cellstr(f);
% filter out unwanted files
for cexcl = 1:length(excluded)
    idx = strfind(cellstr(f),excluded{cexcl});
    for cidx = 1:length(idx)
        if ~isempty(idx{cidx})
            f{cidx} = '';
        end
    end
end

% Table of content
fprintf(fp,'\n### Content\n\n');
for cf = 1:size(f,1)
    if ~isempty(f{cf})
        fnam = deblank(strrep(f{cf},[hMRIdir filesep],''));
        flink = strrep(fnam,' ','-');
        flink = strrep(flink,filesep,'');
        flink = strrep(flink,'.','');
        fprintf(fp,'- [%s](HelpScripts#%s)    \n',fnam,flink);
    end
end
    
% Content
for cf = 1:size(f,1)
    if ~isempty(f{cf})
        fnam = deblank(strrep(f{cf},[hMRIdir filesep],''));
        fprintf(fp,'\n%s\n```\n%s```\n',header(fnam,4),get_mfile_help(f{cf}));
    end
end

fclose(fp);

end

%==========================================================================
function hstr = get_mfile_help(f)
% delimitation must be refined and/or the way to write the help section
% must be more uniform... Should be
% - continuous comment section ("% " at the beginning of each and every
%   line to be included in the Help)
% - ending with copyright or simply an interruption in the "% " block.

hstr = [];

% help text, minus copyrights
htxt = textscan(fopen(f),'%s','delimiter','\n','whitespace','');
htxt = htxt{1};
i_beg = find(strncmp('% ',htxt,2));
i_end = find(strncmp('% Copyright (C)',htxt,15))-2;
if isempty(i_end)
    i_end = i_beg(find(diff(i_beg)>1,1));
end
if isempty(i_beg)
    return
end
i_beg = i_beg(1);
htxt = htxt(i_beg:i_end);

for jj=1:numel(htxt)
    hstr = [hstr sprintf('%s \n',htxt{jj}(2:end))]; %#ok<AGROW>
end

end


%==========================================================================
function ok = chapter(c,fp)

fprintf(fp, '\n%s\n\n', header(c.name,3));
write_help(c,fp);

switch class(c),
    case {'cfg_branch','cfg_exbranch'},
        for i=1:numel(c.val),
            section(c.val{i},fp);
        end;
    case {'cfg_repeat','cfg_choice'},
        for i=1:numel(c.values),
            section(c.values{i},fp);
        end;
end;
ok = true;
end

%==========================================================================
function section(c,fp,lev)

if nargin<3, lev = 1; end;
fprintf(fp,'\n%s- %s    \n',repmat('  ',1,lev-1),bold(c.name));
write_help(c,fp,lev);

switch class(c),
    case {'cfg_branch','cfg_exbranch'},
        for i=1:numel(c.val),
            section(c.val{i},fp,lev+1);
        end
    case {'cfg_repeat','cfg_choice'},
        for i=1:numel(c.values),
            section(c.values{i},fp,lev+1);
        end
end

end

%==========================================================================
function write_help(hlp,fp,lev)

if nargin<3;lev = 1; end

if isa(hlp, 'cfg_item'),
    if ~isempty(hlp.help),
        hlp = hlp.help;
    else
        return;
    end
end

if iscell(hlp),
    for i=1:numel(hlp),
        write_help(hlp{i},fp,lev);
    end
    return
end

str = clean_txt(hlp);
fprintf(fp,'%s%s\n',repmat('  ',1,lev),str);
end

%==========================================================================
function str = clean_txt(str)
str  = strrep(str,'<','');
str  = strrep(str,'>','');
str  = strrep(str,'``','"');
str  = strrep(str,'''','"');
end

%==========================================================================
function str = bold(str)
str = sprintf('**%s**',str);
end

%==========================================================================
function str = italic(str)
str = sprintf('*%s*',str);
end

%==========================================================================
function str = header(str,level)
str = sprintf('%s %s',repmat('#',1,level),str);
end
