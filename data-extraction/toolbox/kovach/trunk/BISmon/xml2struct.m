
function out = xml2struct(xml)

%%% Convert xml to a struct with fields determined by tag names

% C. Kovach 2016

xml = regexprep(xml,'<([^>/\s]*?)\s*/\s*>','<$1></$1>');
parsetext = regexp(xml,'<\s*([^>\s]*)([^>]*)>(.*?)</\1>','tokens');

if isempty(parsetext)
    out = {xml};
else

    parsetext = cat(1,parsetext{:});
    q = unique(parsetext(:,1)');
    q(2,:) = {struct([])};
    out.xmlTagAttributes = struct(q{:});
    for k = 1:size(parsetext,1)
        tagattrib = regexp(parsetext{k,2},'(\w*)="([^"]*)','tokens');
        tagattrib = cat(1,tagattrib{:});
%         out.xmlTagAttributes.(parsetext{k,1})(end+1) = struct;
        for kk =1:size(tagattrib,1)
             out.xmlTagAttributes.(parsetext{k,1})(k).(tagattrib{kk,1}) = tagattrib{kk,2}; 
        end
        if isfield(out,parsetext{k,1}) && isstruct(out.(parsetext{k,1}))
            out.(parsetext{k,1})(end+1) = xml2struct(parsetext{k,3});
        else
           out.(parsetext{k,1}) = xml2struct(parsetext{k,3});
        end   
    end
end



