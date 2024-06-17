function [txt,rvdat,res] = readPage(pagename)

% txt = readPage(pagename)
% Reads wikitext from the page given by pagename. 
%
% See also WRITETOPAGE and MWAPI

% C. Kovach 2017

res = mwapi(struct('action','parse','page',pagename,'prop','wikitext|revid'));

if isfield(res,'parse')
    txt = res.parse.wikitext.x_;
    rvdat = mwapi(struct('action','query','prop','revisions','revids',num2str(res.parse.revid)));
    rvdat = rvdat.query.pages.(sprintf('x%i',res.parse.pageid)).revisions;

elseif ischar(res)
    txt = regexp(res,'"wikitext":{"\*":"(.*?)"}','tokens','once');
    if ~isempty(txt)
        txt = txt{1};
    end
    rvdat='';
else
    txt = '';
    rvdat='';
end