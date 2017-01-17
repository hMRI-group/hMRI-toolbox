function [nFieldFound, fieldList] = findFieldName(inStruct, fieldName, varargin)
% PURPOSE
% To search a structure recursively to retrieve fields with a
% specific name (or part of that name). Returns the number of fields found
% with that name (nFieldFound) and a cell array containing the path of
% each occurrence of that name (each line in fieldList gives the list of
% nested fields down to the expected one, e.g.
% inStruct.(fieldList{1,1}).(fieldList{1,2}). ...
% .(fieldList{1,last-non-empty-one})). If unsuccessful, the 
% search returns [0,{}].  
%
% ARGUMENTS
% - inStruct is the name of the structure to search. If a subfield is a
%   cell array, each element is searched in case it would be structure.
% - fieldName is the name of the field to be found
% - varargin is a series of optional parameters:
%    - matchType: either 'exact' (field name == fieldName) or 'partial'
%      (field name contains fieldName). Default is 'partial'.
%    - caseSens: whether the search is case sensitive or not. caseSens =
%      'sensitive' or 'insensitive'. Default is 'insensitive'.
%
% SYNTAX
% Except for the first two mandatory parameters which must be given in the
% right order, the other arguments (varargin) are passed as pairs of
% parameter name/parameter value, e.g.
% [iFR, fL] = findFieldName(inStruct,'RepetitionTime','matchType','exact')
% [iFR, fL] = findFieldName(inStruct,'Repet','matchType','partial','caseSens','sensitive')
% [iFR, fL] = findFieldName(inStruct,'repet','caseSens','insensitive')
%
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - November 2014 - Cyclotron Research Centre
%--------------------------------------------------------------------------

% sort out the input arguments:
p = inputParser;
p.addRequired('inStruct',@isstruct);
p.addRequired('fieldName',@ischar);
p.addOptional('fieldList',{},@(x) iscell(x));
expectedMatch = {'exact','partial'};
p.addOptional('matchType','partial',@(x) any(validatestring(x,expectedMatch)));
expectedCase = {'insensitive','sensitive'};
p.addOptional('caseSens','insensitive',@(x) any(validatestring(x,expectedCase)));

p.parse(inStruct, fieldName, varargin{:});
   
inStruct = p.Results.inStruct;
fieldName = p.Results.fieldName;
fieldList = p.Results.fieldList;
matchType = p.Results.matchType;
caseSens = p.Results.caseSens;

nFieldFound = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
    switch caseSens
        case 'insensitive'
            s1 = lower(f{i});
            s2 = lower(strtrim(fieldName));
        case 'sensitive'
            s1 = f{i};
            s2 = strtrim(fieldName);
    end
    switch matchType
        case 'exact'
            ok = strcmp(s1,s2);
        case 'partial'
            ok = strfind(s1,s2);
    end
    if ok
        nFieldFound = nFieldFound + 1;
        fieldList{end+1,1} = f{i};
    end
    if iscell(inStruct(1).(f{i}))
        for cc =1:length(inStruct(1).(f{i}))
            if isstruct(inStruct(1).(f{i}){cc})
                [iFR, fL] = findFieldName(inStruct(1).(f{i}){cc}, fieldName, ...
                    'matchType', matchType, 'caseSens', caseSens);
                if iFR
                    nFieldFound = nFieldFound + iFR;
                    for ci = 1:size(fL,1)
                        fieldList{end+1,1} = sprintf('%s{%d}',f{i},cc);
                        for cj = 1:size(fL,2)
                            fieldList{end,1+cj} = fL{ci,cj};
                        end
                    end
                end
            end
        end
    elseif isstruct(inStruct(1).(f{i}))
        [iFR, fL] = findFieldName(inStruct(1).(f{i}), fieldName, ...
            'matchType', matchType, 'caseSens', caseSens);
        if iFR
            nFieldFound = nFieldFound + iFR;
            for ci = 1:size(fL,1)
                fieldList{end+1,1} = f{i};
                for cj = 1:size(fL,2)
                    fieldList{end,1+cj} = fL{ci,cj};
                end
            end
        end
    end
end