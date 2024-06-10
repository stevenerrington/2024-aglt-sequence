function res = writeToPage(me,page_name,text,action,comment,tag)

%
% Usage:
%
%   mwapi.result = writeToPage(page_name,text)
%
% Appends text to page with title given by page_name
%
%
%   mwapi.result = writeToPage(page_name,text,action)
%
% Append, prepend or overwrite, depending on the value of action:
%   action = -1 :  prepend
%   action = 0  :  overwrite
%   action = 1  :  append (default)
%
%See also READPAGE and MWAPI

% C. Kovach 2017

if nargin < 4 || isempty(action)
    action = 1;
end
if nargin <6 || isempty(tag)
    tag='lwapi-edit';
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


text = regexprep(text,'&','%26');
reqstr.(fld) = text;
if nargin >4 && ~isempty(comment)
    reqstr.summary = comment;
end

tok = me.get_token();
reqstr.token = tok;
reqstr.format='json';

res = me.request(reqstr);
if ischar(res)
    res = me.fromJSON(res);
end
if isfield(res,'edit') && strcmpi(res.edit.result,'Success') && ~isfield(res.edit,'nochange')
    try
    tagstr = struct('action','tag','revid',num2str(res.edit.newrevid),'add',tag,'token',tok,'format','json');
    res2 = me.request(tagstr);
    catch 
    end
end
