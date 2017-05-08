function hmri_latex(c)
% Extract information from the toolbox m-files and output them as usable
% .tex files which can be directly included in the manual.
%
% There are 2 types of m2tex operations:
% 1. converting the job configuration tree, i.e. *_cfg_* files defining the
%    batching interface into a series of .tex files.
% 2. converting the help header of the functions into .tex files.
%
% These files are then included in a manually edited hmri_manual.tex file,
% which also includes chapter/sections written manually.
%
% File derived from that of the SPM8 distribution.
% http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by John Ashburner & Christophe Phillips

skip_fct = true;
% If true, then the list of function is not processed and the corresponding
% .tex file not created.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Turning the cfg files into a .tex file
if ~nargin,
    c = tbx_cfg_hmri;
end

hMRIdir = spm_str_manip(mfilename('fullpath'),'h');

for ii=1:numel(c.values),
    bn = c.values{ii}.tag;
    fp = fopen(fullfile(hMRIdir,'manual',['batch_',bn,'.tex']),'w');
    if fp==-1, sts = false; return; end;
    chapter(c.values{ii},fp);
end;

if ~skip_fct
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. picking all the functions help files and put them into functions.tex
    fp = fopen(fullfile(hMRIdir,'manual','adv_functions.tex'),'w');
    if fp==-1, sts = false; return; end;
    l_subdirs = {};
    excl_files = {};
    
    % Heading part
    fprintf(fp,'\\chapter{%s  \\label{Chap:%s}}\n\\minitoc\n\n\\vskip 1.5cm\n\n',...
        texify('List of hMRI functions'),'sec:functions');
    fprintf(fp,'This is the list of hMRI functions, including the subdirectories: ');
    for ii=1:numel(l_subdirs)
        fprintf(fp,'%s',texify(['{\tt ',l_subdirs{ii},'}']));
        if ii<numel(l_subdirs)-1
            fprintf(fp,', ');
        elseif ii==numel(l_subdirs)-1
            fprintf(fp,' and ');
        else
            fprintf(fp,'.\n\n');
        end
    end
    
    % Deal with main directoy 1st
    f = spm_select('List',PRTdir,'.*\.m$');
    for ii=1:numel(excl_files)
        f(strcmp(cellstr(f),excl_files{ii}),:) = [];
    end
    write_mfiles_help(f,fp);
    
    % Deal with subdirectories
    for ii = 1:numel(l_subdirs)
        p = fullfile(PRTdir,l_subdirs{ii});
        fprintf(fp,'\n\\%s{%s}\n','section',texify(l_subdirs{ii}));
        f = spm_select('List',p,'.*\.m$');
        write_mfiles_help(f,fp,l_subdirs(ii));
    end
end
end

%==========================================================================
function write_mfiles_help(f,fp,base_dir)

if nargin<3,
    base_dir = '';
    lev = 1;
else
    lev = numel(base_dir)+1;
    tmp = base_dir{1};
    ltmp = tmp;
    for ii=2:numel(base_dir)
        ltmp = [tmp,'\textbackslash ',base_dir{ii}];
        tmp = fullfile(tmp,base_dir{ii});
    end
    lbase_dir = ltmp;
    base_dir = tmp;
end

sec = {'section','subsection','subsubsection','paragraph','subparagraph', ...
    'textbf','textsc','textsl','textit'};

for ii=1:size(f,1)
    % section
    if isempty(base_dir)
        lfunc_name = deblank(f(ii,:));
    else
        lfunc_name = [lbase_dir,'\textbackslash ',deblank(f(ii,:))];
    end
    func_name = fullfile(base_dir,deblank(f(ii,:)));
    fprintf(fp,'\n\\%s{%s}\n',sec{min(lev,length(sec))},texify(lfunc_name));
    fprintf(fp,'%s\n\n',texify('\begin{alltt}'));
    
    % help text, minus copyrights
    htxt = textscan(fopen(func_name),'%s','delimiter','\n','whitespace','');
    htxt = htxt{1};
    i_beg = find(strncmp('% ',htxt,2)); i_beg = i_beg(1);
    i_end = find(strncmp('% Copyright (C)',htxt,15))-2;
    htxt = htxt(i_beg:i_end);
    
    for jj=1:numel(htxt)
        if strcmp(htxt{jj},'%')
            fprintf(fp,'%s\n',' ');
        else
            fprintf(fp,'%s\n',texify(htxt{jj}(2:end)));
        end
    end
    fprintf(fp,'%s\n\n',texify('\end{alltt}'));
end

end

%==========================================================================
function sts = chapter(c,fp)
bn = c.tag;
if nargin<2
    fp = fopen(fullfile(pwd,'manual',[bn,'.tex']),'w');
    if fp==-1, sts = false; return; end;
end

fprintf(fp,'%% $Id: hmri_latex.m 31 2013-11-27 16:42:58Z christophe $ \n\n');
fprintf(fp, ...
    '\\chapter{%s  \\label{Chap:%s}}\n\n\\vskip 1.5cm\n\n', ...
    texify(c.name),c.tag);
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
fclose(fp);
sts = true;
end

%==========================================================================
function section(c,fp,lev)
if nargin<3, lev = 1; end;
sec = {'section','subsection','subsubsection','paragraph','subparagraph', ...
    'textbf','textsc','textsl','textit'};
% if lev<=length(sec),
fprintf(fp,'\n\\%s{%s}\n',sec{min(lev,length(sec))},texify(c.name));
write_help(c,fp);
switch class(c),
    case {'cfg_branch','cfg_exbranch'},
        for i=1:numel(c.val),
            section(c.val{i},fp,lev+1);
        end;
    case {'cfg_repeat','cfg_choice'},
        for i=1:numel(c.values),
            section(c.values{i},fp,lev+1);
        end;
end;
% else
if lev>length(sec),
    warning(['Too many nested levels... ',c.name]); %#ok<WNTAG>
end
end

%==========================================================================
function write_help(hlp,fp)
if isa(hlp, 'cfg_item'),
    if ~isempty(hlp.help),
        hlp = hlp.help;
    else
        return;
    end;
end;
if iscell(hlp),
    for i=1:numel(hlp),
        write_help(hlp{i},fp);
    end;
    return;
end;
str = texify(hlp);
fprintf(fp,'%s\n\n',str);
end

%==========================================================================
function str = texify(str0)
st1  = strfind(str0,'/*');
en1  = strfind(str0,'*/');
st = [];
en = [];
for i=1:numel(st1),
    en1  = en1(en1>st1(i));
    if ~isempty(en1),
        st  = [st st1(i)];
        en  = [en en1(1)];
        en1 = en1(2:end);
    end;
end;

str = [];
pen = 1;
for i=1:numel(st),
    str = [str clean_latex(str0(pen:st(i)-1)) str0(st(i)+2:en(i)-1)];
    pen = en(i)+2;
end;
str = [str clean_latex(str0(pen:numel(str0)))];
end

%==========================================================================
function str = clean_latex(str)
str  = strrep(str,'$','\$');
str  = strrep(str,'&','\&');
str  = strrep(str,'^','\^');
str  = strrep(str,'_','\_');
str  = strrep(str,'#','\#');
% str  = strrep(str,'\','$\\$');
str  = strrep(str,'|','$|$');
str  = strrep(str,'>','$>$');
str  = strrep(str,'<','$<$');
end

%==========================================================================
function bibcstr = get_bib(bibdir)
biblist = dir(fullfile(bibdir,'*.bib'));
bibcstr={};
for k = 1:numel(biblist)
    [p n e v] = spm_fileparts(biblist(k).name);
    bibcstr{k}  = fullfile(bibdir,n);
end
end

%==========================================================================
function clean_latex_compile
PRTdir = prt('dir');
p = fullfile(PRTdir,'manual');
[f, d] = spm_select('FPlist',p,'.*\.aux$');
f = strvcat(f, spm_select('FPlist',p,'.*\.tex$'));
f = strvcat(f, spm_select('FPlist',p,'^manual\..*$'));
f(strcmp(cellstr(f),fullfile(PRTdir,'manual','prt_manual.tex')),:) = [];
f(strcmp(cellstr(f),fullfile(PRTdir,'manual','prt_manual.pdf')),:) = [];
for i=1:size(d,1)
    f = strvcat(f, spm_select('FPlist',deblank(d(i,:)),'.*\.aux$'));
end
f(strcmp(cellstr(f),filesep),:) = [];
disp(f); pause
for i=1:size(f,1)
    spm_unlink(deblank(f(i,:)));
end
end
