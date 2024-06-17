
function [out,graph,nfp] = dependencytree(in,get_file_content,maxdepth,out,depth,nfp)

% [out,graph] = dependencytree(mfile)
%
% Creates a directed graph of dependencies for a given file or set of files.
%
% Inputs:
%   mfile - filename or cell array of filenames to analyze.
%
% Outputs:
%   out - struct array containing dependencies of mfile with fields:
%      .filename - name of the file
%      .dependencies - Files that out.filename depends on, as indices into out
%      .serves  -  All files in out that depend on out.filename.
%
%   graph - A matlab digraph object with the graph structure. 
%


if ischar(in)
    in = struct('filename',which(in),'dependson',[],'serves',[],'index',[],'content','');
elseif iscell(in)
    in = struct('filename',in,'dependson',[],'serves',[],'index',[],'content','');
% else
%     out = in(1:end-1);
%     in = out(end);
end

if nargin < 3 || isempty(maxdepth)
   maxdepth = Inf; 
end

if nargin < 2 || isempty(get_file_content)
    get_file_content = true;
end

if nargin < 5 || isempty(depth)
    depth = 0;
end

if nargin < 6 || isempty(nfp)
    nfp = 0;
end

if nargin < 4 || isempty(out)
    out = struct('filename',{},'dependson',{},'serves',{},'index',{},'content',{});
end

for k = 1:length(in)
    fns = {out.filename};
    in(k).filename = which(in(k).filename);
    if get_file_content
        in.content = '';
    end
    if ismember(in(k).filename,fns)
        continue
    else 
        out(end+1) = in(k);
        indx = length(out);
        out(end).index = indx;
        if get_file_content
            fid = fopen(out(end).filename,'r');
            out(end).content = fread(fid,'uchar=>char')';
            fclose(fid);
        end
        nfp = fprintf([repmat('\b',1,nfp),'\nDepth %i: scanning %s ...'], depth,in(k).filename)-nfp;

        if maxdepth > 1
            deps = matlab.codetools.requiredFilesAndProducts(in(k).filename,'toponly');
        else
            deps = matlab.codetools.requiredFilesAndProducts(in(k).filename);
        end
        [~,out(end).dependson] = ismember(deps,fns);
        missing = find(~out(end).dependson);
        
        for kk = 1:length(missing)
            [a,b] = ismember(deps{missing(kk)},{out.filename});            
            if a
                out(indx).dependson(missing(kk)) = b;
            elseif maxdepth <= 1
                out(end+1).filename = dep{missing(kk)};
                depinx = length(out);
                out(end).dependson = depindx;
                out(end).index = depindx;
                out(indx).dependson(kk) = depindx;
              
            else
                depindx = length(out);
                [out,~,nfp] = dependencytree(deps{missing(kk)},get_file_content,maxdepth-1,out,depth+1,nfp);
                out(indx).dependson(missing(kk)) = depindx+1;
            end
        end
    end    
end

if depth ==0
    for k = 1:length(out)
        out(k).serves = find(arrayfun(@(x)any(x.dependson==k),out));
    end
    if nargout > 1
        A = arrayfun(@(x)cat(1,x.serves,ones(size(x.serves))*x.index),out,'uniformoutput',false); A = [A{:}];
        B = arrayfun(@(x)cat(1,ones(size(x.dependson))*x.index,x.dependson),out,'uniformoutput',false); B= [B{:}];
        AB = unique(cat(2,A,B)','rows');
        
        [~,fn] = fileparts({out.filename});
        fn = regexprep(fn,'_','\\_');
        graph = digraph(fn(AB(:,1)),fn(AB(:,2)),'omitselfloops');
        
    end
else
    graph = [];
end
