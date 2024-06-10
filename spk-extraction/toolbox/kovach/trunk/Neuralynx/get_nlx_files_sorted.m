function [fnames,fsizes] = get_nlx_files_sorted(get_path,pth,chan_types)


%%% Return nlx and nev files sorted numerically 

if nargin < 1  ||  isempty(get_path)
    get_path = false;
end

if nargin < 3 || isempty(chan_types)
    chan_types = {''};
else
    if ischar(chan_types)
        chan_types = {chan_types};
    end
end


persistent path 

if nargin > 1 && ~isempty(pth)
    path = pth;
end

fnames = {};
fsizes=[];
for k = 1:length(chan_types)
    if ~get_path 
        [fns,pth] = uigetfile({strcat(chan_types{k},'*.ncs;*.nev'),'NLX files (*.ncs,*.nev)';'*.nev','Event files only'},'Select NLX data file',path,'multiselect','on');
        if isnumeric(fns)
            return
        end
        if ~iscell(fns)
            fns = {fns};
        end
    else
        if nargin < 2 && isempty(pth)
            pth = uigetdir(path,'Select NLX data file directory');
        end
        fns = cat(1,dir(fullfile(pth,[chan_types{k},'*.ncs'])),dir(fullfile(pth,'*.nev')));
        fsz=[fns.bytes];
        fns = {fns.name};
    end
    if isnumeric(fns)
         return
    else
        path = pth;
    end
  



    % sort numerically
    decell= @(x)[x{:}];
    reg = regexp(fns,'(\d+)_*(\d*)','tokens','once');
    reg(cellfun(@isempty,reg)) = {{'',''}};
    chtypes = regexp(fns,'([^\d]+)','tokens','once');
    chtypes(cellfun(@isempty,chtypes)) = {{''}};
    chtypes =cat(1,chtypes{:});
    [~,~,unqi] = unique(chtypes);
    
    chn = fliplr(arrayfun(@(x)str2double(['0',decell(x)]),cat(1,reg{:})));
    [~,srti] = sortrows([unqi,chn]);

    fns = fullfile(path,fns(srti));

    if get_path
        fsz = fsz(srti);
    else
        fsz = nan(size(fns));
    end

    if ~iscell(fns)
        fns = {fns};
    end

    fnames = [fnames,fns]; %#ok<*AGROW>
    fsizes = [fsizes,fsz]; %#ok<*AGROW>
end
[fnames,unqi] = unique(fnames,'stable');
fsizes=fsizes(unqi);
