
rmpath('~/svn/KovachToolbox/GazeReader/0.1/')
 dat = readtable('allblocks2.csv','delimiter',',');
% mw = mwapi2;
% addpath('/home/kovachc/svn/KovachToolbox/trunk/GazeReader/0.1/')
dat.HasSubprotocol = upper(dat.HasSubprotocol);
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Tone Mapping','IgnoreCase',true)) = {'Tone Mapping'};
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Tone Noise Target','IgnoreCase',true)) = {'Tone Noise Target'};
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Stream2','IgnoreCase',true)) = {'Stream2'};
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Basic Click Map','IgnoreCase',true)) = {'Basic Click Mapping'};
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Speech Noise','IgnoreCase',true)) = {'Speech Noise'};
dat.HasSubprotocol(contains(dat.HasSubprotocol,'Random Click Rates','IgnoreCase',true)) = {'Random Click Rates'};

[protos,~,proti] = unique(dat.HasSubprotocol);
[subs,~,subi] = unique(dat.HasSubject);

ct = crosstab(subi,proti);
%Choose the largest cluster of subjects with common protocols 
CC = corr((ct>0)');
[scclass,e,v] = spectralcluster(CC,5);
[srt,srti] = sort(sum(ct>0),'descend');
scnum = hist(scclass,1:10);
[mx,getsubi] = max(scnum);
getsub = find(scclass==getsubi);
CC = (ct(getsub,:)>0)'*(ct(getsub,:)>0);
CC(isnan(CC))=0;
[u,l,v] = svd(CC);
[srt,srti] = sort(u(:,1)*sign(mean(u(:,1))),'descend');

getprotos = [srti(1:50); find(contains(protos,'value attention','ignorecase',true)); find(contains(protos,'morphlet','ignorecase',true))];

sf = @(x) ~cellfun(@isempty,strfind(lower(protos),lower(x)));
dropprot = sf('abort') | sf('fMRI')| sf('EsTT')| sf('Es-Amyg')| sf('EsAmyg')|sf('baseline')|cellfun(@isempty,protos);

getbl =find(ismember(subi,getsub)&ismember(proti,getprotos)&~dropprot(proti));
% reseed(12345)
% reseed(54321)
reseed(45312)
getbl = getbl(randperm(length(getbl)));

%%% Added after pseudorandom shuffling
dropprot2 = sf('sleep')|sf('ESTT')|sf('ES TT')|(sf('STIM') & (sf('ES')| sf('ET')));
% dropbl = ~ismember(subs(subi(getbl)),{'Subject/247','Subject/422'});
getbl(dropprot2(proti(getbl))) = [];

getbl(ismember(dat.Var1(getbl),{'250-060','287-050','232-070'}))=[]; %These are stimulation blocks
% Matrix of subjects x protocls 

% Limit blocks of the same protocol to 2 per subject
maxn = 2;
[unq,~,unqi] = unique([subi(getbl),proti(getbl)],'rows');
blmat = zeros(length(getbl),length(unq));
blmat((1:length(getbl))' + (unqi-1)*length(getbl))=1;
blmat = blmat.*(cumsum(blmat)<=maxn);
getbl = getbl(sum(blmat,2)>0);

run_blocks= dat.Var1(getbl);

% run_blocks=run_blocks(randperm(length(run_blocks)));
if ~exist('submit_profiles','var')
submit_profiles = struct('profile',{'mid_mem','high_mem','std_mem'},...
                         'cluster',{'argon','argon','argon'},...
                         'nslots',{4,2,8},...
                         'queue',{'CCOM,UI,all.q','CCOM,UI-HM,NEUROSURGERY,all.q','CCOM,UI,all.q'},...
                         'xnes',{{}});
end

opts.ncomp = 3;
opts.windur = sqrt(17);
opts.do_regression = false;
opts.lowpass=200;
opts.glowpass=200;
opts.save_space = true;
opts.run_phase_randomized=true;
opts.run_compiled = false;
local_data_dir = 'the_big_analysis';
mdl =[];
%%
% for repeat = 1:10
    for k = 1:length(run_blocks)
        if exist('blkss','var') && length(blkss)>=k && strcmp(blkss(k).block,run_blocks{k})
            blks = blkss(k);
        else
            blks = get_protocol_labwiki(run_blocks{k});
            try
                blks = locateNlx(blks);
            catch
            end
            if isfield(blks,'blkfiles') && ~isempty(blks.blkfiles)
                            blkss(k) = blks;
            else
                fprintf('\nFailed to find data for %s',run_blocks{k});
                continue
            end
        end
        try
        bispectral_analysis
        catch err
            warning(err.message)
        end
        
        xness{k} = xnes;
        k
    end
    fprintf('\nPausing for a bit...')
    pause(1000)
% end

%%

for repeat = 1:1000
    try
    jb = get_jobs;
    q = jb(arrayfun(@(x)contains(x.local_save_dir,'the_big')&&isequal(x.status,'finished'),jb));
    
%     q = jb(arrayfun(@(x)isequal(x.status,'error') && contains(x.local_save_dir,'the_big'),jb));
%     q = q(str2double({q.jobid})<1407764);
   %%
    re = regexp({q.local_save_dir},'the_big_analysis/(\d\d\d-\d\d\d)','tokens','once'); 
    run_blocks2=[re{:}]
    for k = 1:length(run_blocks2)
        
        blks = get_protocol_labwiki(run_blocks2{k}) ;
        try
         bispectral_analysis
         if already_done
             q(k).close();
         end
        catch err
            warning(err.message)
        end
       k
    end       
    catch err
        warning(err.message)
    end
   fprintf('\nPausing for a bit...')
    pause(60)
end

%%
d = dir(local_data_dir);
re = regexp({d.name},'\d\d\d-\d\d\d','match','once');
gotblock = find(ismember(dat.Var1,re));
[spmat,~,~,lbl] = crosstab(subs(subi(gotblock)),protos(proti(gotblock)));
figure, 
sz= size(spmat);
plot([[0 sz(1)+1]'*ones(1,sz(2)+1),[1 1]'*(.5:sz(1)+.5)],[[1 1]'*(.5:sz(2)+.5),[0 sz(2)+1]'*ones(1,sz(1)+1)],'color',[1 1 1]/2)
[asrt,asrti] = sort(lbl(1:size(spmat,1),1));
[bsrt,bsrti] = sort(lbl(1:size(spmat,2),2));
[a,b] = ind2sub(sz,find(spmat(asrti,bsrti)));
hold on, plot(a,b,'k*')
set(gca,'ytick',1:sz(2),'yticklabel',lbl(bsrti,2),'xtick',1:sz(1),'xticklabel',regexprep(lbl(asrti,1),'Subject/',''),'xticklabelrotation',90)
axis([-1 sz(1)+.5 -1 sz(2)+.5])
axis image

[srt,srti] = sort(sum(spmat>0,2));

getbl2 = getbl(ismember(subs(subi(getbl)),asrt(srti(srt<3))));






