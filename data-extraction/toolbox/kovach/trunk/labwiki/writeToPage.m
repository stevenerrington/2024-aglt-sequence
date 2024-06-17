function res = writeToPage(page_name,text,action,comment)

%
% Usage:
%
%   result = writeToPage(page_name,text)
%
% Appends text to page with title given by page_name
%
%
%   result = writeToPage(page_name,text,action)
%
% Append, prepend or overwrite, depending on the value of action:
%   action = -1 :  prepend
%   action = 0  :  overwrite
%   action = 1  :  append (default)
%
%See also READPAGE and MWAPI

% C. Kovach 2017

if nargin < 3 || isempty(action)
    action = 1;
end

switch action
    case {-1,'prepend'}
        fld = 'prependtext';
    case {0,'overwrite'}
        fld = 'text';
    case {1,'append'}
        fld = 'appendtext';
    otherwise
        error('Unrecognized action. Must be prepend (-1), overwirte (0), or append (1).')
end

reqstr.title = page_name;
reqstr.action = 'edit';

reqstr.(fld) = text;
if nargin >3 && ~isempty(comment)
    reqstr.summary = comment;
end
    
res = mwapi(reqstr);
if isfield(res,'edit') && strcmpi(res.edit.result,'Success')
    tagstr = struct('action','tag','revid',num2str(res.edit.newrevid),'add','lwapi-edit');
    res2 = mwapi(tagstr);
end
