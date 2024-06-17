function [txt,rvdat,res] = readPage(me,pagename)

% txt = readPage(pagename)
% Reads wikitext from the page given by pagename. 
%
% See also WRITETOPAGE and MWAPI

% C. Kovach 2017

res = me.request(struct('action','parse','page',pagename,'prop','wikitext|revid','format','json'));

if ischar(res)
    res = me.fromJSON(res);
end

if isfield(res,'parse')
    txt = res.parse.wikitext.x_;
    if nargout >1
        rvdat = me.request(struct('action','query','prop','revisions','revids',num2str(res.parse.revid),'format','json'));
        if ischar(rvdat)
            rvdat = me.fromJSON(rvdat);
        end
        rvdat = rvdat.query.pages.(sprintf('x%i',res.parse.pageid)).revisions;
    end
else
    txt = '';
    rvdat='';
end