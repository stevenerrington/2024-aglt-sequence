function xml = struct2xml(str)

% Converts a structure whose fields contain either substructures or
% character arrays into xml text. The result is xml text with tags given 
% by the structure field names and attributes contained in the special field,
% 'xmlTagAttributes'.
%
% See also XML2STRUCT

% C. Kovach 2016

fld = setdiff(fieldnames(str),'xmlTagAttributes');
xml = '';
for k = 1:length(fld)
    if ~isstruct(str.(fld{k})) && ~iscell(str.(fld{k}))
      str.(fld{k}) = {str.(fld{k})};
    end
    for j = 1:length(str.(fld{k}))
        if isfield(str,'xmlTagAttributes') && isfield(str.xmlTagAttributes,fld{k})
            tagparams = fieldnames(str.xmlTagAttributes.(fld{k})(j))';
            tagparams(2,:) = cellfun(@(x)str.xmlTagAttributes.(fld{k})(j).(x),tagparams,'uniformoutput',false);
            tagstr = deblank(sprintf(' %s="%s"',tagparams{:}));
        else
            tagstr = '';
        end

        head = sprintf('\t\n<%s%s>',fld{k},tagstr);
        foot = sprintf('</%s>',fld{k});
        if isstruct(str.(fld{k}))
            body = regexprep(struct2xml(str.(fld{k})(j)),'\n','\n\t');
        else
            body = str.(fld{k}){j};
        end
        xml = sprintf('%s%s%s\n%s',xml,head,body,foot);
    end
end
