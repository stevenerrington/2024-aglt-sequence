
function xne = xnemomdist(X,Y,order)

if isscalar(Y)
    order = Y;
    Y = X;
end

xne = xargon('momdist_clust');

out.order = order;
out.X = X;
out.Y=Y;
out.chunk = 5;

tempfn = fullfile(xne.tempdir,'input.mat');
save(tempfn,'-struct','out');

xne.datafiles = tempfn;

xne.queue = 'all.q,UI,CCOM';
xne.profile = 'mid_mem';
xne.nslots = 2;

xne.nparallel = ceil(size(Y,2)./out.chunk);


xne.finish = @(varargin)fin(xne);

xne.run_compiled=true;
xne.create_path();
xne.compile();
pause(1);
xne.submit();

function fin(xne)

switch  xne.status
    
    case 'finished'
        d = dir(fullfile(xne.subpaths.output.local,'out*.mat'));
        re = regexp({d.name},'out(\d*)[.]mat','tokens','once');
        re = cellfun(@str2double,re);
        xne.jobindices = setdiff(xne.jobindices,re);
        if isempty(xne.jobindices)
            xne.default_finish();
            d = dir(fullfile(xne.local_save_dir,'out*.mat'));
            re = regexp({d.name},'out(\d*)[.]mat','tokens','once');
            [re,srti] = sort(cellfun(@str2double,re));
            d = d(srti);
            for k = 1:length(re)
                outs(k) = load(fullfile(d(k).folder,d(k).name));
            end
            out.Dst = cat(1,outs.Dst);
            out.Dsgn = cat(1,outs.Dsgn);
            out.pkcorr = cat(1,outs.pkcorr);
            save(fullfile(xne.local_save_dir,sprintf('result%i.mat',outs(1).order)),'-struct','out');
            delete(fullfile(xne.local_save_dir,'out*.mat'))
        else
            fprintf('\n%i of %i finished',length(re),length(xne.jobindices) + length(re));
            xne.make_bash_script;
            xne.compile
            xne.submit(true);
        end
    otherwise
        fprintf('Status is %s.',xne.status);
end
             


