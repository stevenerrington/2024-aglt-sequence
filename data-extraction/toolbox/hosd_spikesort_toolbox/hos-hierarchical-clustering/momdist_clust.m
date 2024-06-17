function momdist_clust(inputdir,outputdir,inputfiles,jobindex)

savefname = fullfile(outputdir,sprintf('out%i',jobindex));
if exist(savefname,'file')
    return
end


ld = load(fullfile(inputdir,inputfiles{1}));

if ~isfield(ld,'chunk')
    ld.chunk = 1;
end

geti = ld.chunk*(jobindex-1)+1:min(ld.chunk*jobindex,size(ld.X,2));
[out.Dst,out.Dsgn,out.pkcorr] = momdist(ld.X(:,geti),ld.Y,ld.order);
out.order = ld.order;

save(savefname,'-struct','out');

