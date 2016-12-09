function varargout = hmri_jsonwrite(varargin)
% Serialize a JSON (JavaScript Object Notation) structure
% FORMAT spm_jsonwrite(filename,json)
% filename - JSON filename
% json     - JSON structure
%
% FORMAT S = spm_jsonwrite(json)
% json     - JSON structure
% S        - serialized JSON structure (string)
%
% FORMAT [...] = spm_jsonwrite(...,opts)
% opts     - structure of optional parameters:
%              compact: compact vs pretty-print formatting [true]
% 
% References:
%   JSON Standard: http://www.json.org/
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_jsonwrite.m 6947 2016-11-23 16:12:12Z guillaume $

%=========================================================================%
% Version slightly modified for better formatting and readability of the
% JSON metadata (when opts.compact = false):
% - tabs used instead of spaces when opt.compact = false
% - arrays containing numerical values are kept on a single line (rather
%   than fully developed with new line for each element in the array)
%=========================================================================%
% Adapted by Evelyne Balteau - Cyclotron Research Centre - December 2016
%=========================================================================%

%-Input parameters
%--------------------------------------------------------------------------
opts         = struct('compact',true);
opt          = struct([]);
if nargin > 1
    if ischar(varargin{1})
        filename = varargin{1};
        json     = varargin{2};
        root     = inputname(2);
    else
        filename = '';
        json = varargin{1};
        opt  = varargin{2};
        root = inputname(1);
    end
    if nargin > 2
        opt  = varargin{3};
    end
else
    filename = '';
    json     = varargin{1};
    root     = inputname(1);
end
fn = fieldnames(opt);
for i=1:numel(fn)
    opts.(fn{i}) = opt.(fn{i});
end

%-JSON serialization
%--------------------------------------------------------------------------
if ~isstruct(json) && ~iscell(json) && ~isa(json,'containers.Map')
    if ~isempty(root)
        json = struct(root,json);
    else
        error('Invalid JSON structure.');
    end
end
if opts.compact, tab = NaN; else tab = 0; end
S = jsonwrite_var(json,tab);

%-Output
%--------------------------------------------------------------------------
if isempty(filename)
    varargout = { S };
else
    fid = fopen(filename,'wt');
    if fid == -1
        error('Unable to open file "%s" for writing.',filename);
    end
    fprintf(fid,'%s',S);
    fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = jsonwrite_var(json,tab)
if nargin < 2, tab = 0; end
if isstruct(json) || isa(json,'containers.Map')
    S = jsonwrite_struct(json,tab);
elseif iscell(json)
    S = jsonwrite_cell(json,tab);
elseif ischar(json)
    S = jsonwrite_char(json);
elseif isnumeric(json) || islogical(json)
    S = jsonwrite_numeric(json);
else
    error('Class "%s" is not supported.',class(json));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = jsonwrite_struct(json,tab)
if numel(json) == 1
    if isstruct(json), fn = fieldnames(json); else fn = keys(json); end
    S = ['{' fmt('\n',tab)];
    for i=1:numel(fn)
        if isstruct(json), val = json.(fn{i}); else val = json(fn{i}); end
        % hmri changes: tabs instead of spaces when ~isnan(tab)
        if ~isnan(tab)
            S = [S repmat(fmt('\t',0),1,tab+1) jsonwrite_char(fn{i}) ':' fmt(~isnan(tab)) ...
                jsonwrite_var(val,tab+1)];
        else
            S = [S fmt((tab+1)*2) jsonwrite_char(fn{i}) ':' fmt(~isnan(tab)) ...
                jsonwrite_var(val,tab+1)];
        end
        if i ~= numel(fn), S = [S ',']; end
        S = [S fmt('\n',tab)];
    end
    % hmri changes: tabs instead of spaces when ~isnan(tab)
    if ~isnan(tab)
        S = [S repmat(fmt('\t',0),1,tab) '}'];
    else
        S = [S fmt(2*tab) '}'];
    end
else
    S = jsonwrite_cell(arrayfun(@(x) {x},json),tab);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = jsonwrite_cell(json,tab)
for i=1:numel(json)
    % hmri changes: tabs instead of spaces and special case for numerical
    % values so arrays are not fully developped (better readability):
    if isnumeric(json{i}) || islogical(json{i})
        if i==1, S = '[';end
        S = [S jsonwrite_var(json{i},tab+1)];
        if i~=numel(json), S = [S ',']; end
        if i==numel(json), S = [S ']']; end
    else
        if i==1, S = ['[' fmt('\n',tab)];end
        S = [S repmat(fmt('\t',0),1,tab+1) jsonwrite_var(json{i},tab+1)];
        if i~=numel(json), S = [S ',' fmt('\n',tab)]; end
        if i==numel(json), S = [S fmt('\n',tab) repmat(fmt('\t',0),1,tab) ']']; end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = jsonwrite_char(json)
% any-Unicode-character-except-"-or-\-or-control-character
% \" \\ \/ \b \f \n \r \t \u four-hex-digits
json = strrep(json,'\','\\');
json = strrep(json,'"','\"');
S = ['"' json '"'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = jsonwrite_numeric(json)
if numel(json) > 1
    idx = find(size(json)~=1);
    if numel(idx) == 1 % vector
        S = jsonwrite_cell(num2cell(json),NaN);
    else % array
        S = jsonwrite_cell(num2cell(json,setdiff(1:ndims(json),idx(1))),NaN);
    end
    return;
end
if islogical(json)
    if json, S = 'true'; else S = 'false'; end
else
    S = num2str(json);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = fmt(varargin)
b = '';
if nargin == 1
    if ~isnan(varargin{1}), b = blanks(varargin{1}); end
elseif nargin == 2
    if ~isnan(varargin{2}), b = sprintf(varargin{1}); end
end
