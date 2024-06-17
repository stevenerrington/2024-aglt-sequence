function out = parsejson(json)

open = json=='{';
close = [' ',json(1:end-1)]=='}';
boundaries = open-close;
nest = cumsum(boundaries);
sections = cumsum(open.*(nest==1));

for k = 1:max(sections);

    subsections = cumsum(open(sections==k).*(nest(sections==k)==2));
    subjson = json(sections==k);
    subnest = nest(sections==k);
    
  %  data = subjson(subnest==1);
    subjson = regexp(subjson,'{(.*)}','tokens','once');
    subjson=subjson{1};
    re = regexp(subjson,'"\w*?":');
    re(end+1)=length(subjson);
    supern =  subsections([re(2:end) ])-subsections(re(1:end-1)) ;
    for kk = 1:length(re)-1
        txt = subjson(re(kk):re(kk+1));
        re2 = regexp(txt(1:end-1),'"(\w*)":(.*)','tokens','once');
        if supern(kk)>0
            sub = re2{2};
            sub = regexprep(sub,'\\"','"');
            subout = parsejson(sub);
            if strcmp(re2{1},'form_data')
                c = squeeze(struct2cell(subout));
                c(1,:) = regexprep(c(1,:),'\s*','_');
                c(1,:) = regexprep(c(1,:),',','');
                subout = struct(c{:});
            end
        else
            subout = regexprep(re2{2},'"','');
        end
      out(k).(re2{1})=subout; %#ok<AGROW>
    end
end
        
    