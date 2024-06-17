%  inputdir = 'global_feature_anls_clust_out';
inputdir = xne.subpaths.output.local;
d = dir(fullfile(inputdir,'allfeats*.mat'));

clear ld
for k = 1:length(d)
    ld(k) = load(fullfile(inputdir,d(k).name));
    k
end

allfeats = [ld.allfeats];
allfeatsphaserand = [ld.allfeatsphaserand];

sk = [allfeats.skewness];
skpr = [allfeatsphaserand.skewness];
[p,b] = hist(sk(:),-1:.001:1);
[ppr,b] = hist(skpr(:),-1:.001:1);
thr = b(find(diff(cumsum(ppr)/sum(ppr)<.95)));

figure, plot(sort(sk(:),'descend'))
hold on, plot(sort(skpr(:),'descend'))
legend({'Data','Phase randomized data'})
xlabel 'feature rank'
ylabel 'Filtered data skewness'
set(gca,'yscale','log','xscale','log')
grid on
ylim([.01 100])
axis tight

%%
chunk = 1e3;
fk = 1:chunk:length(selectedfeats);
mkdir('selectedfeats')
clear out
for k = fk
    inds = k-1+(1:chunk);
    out.chunk = selectedfeats(inds);
    out.chunkphrand = selectedfeatsphrand(k-1+(1:chunk));
    out.inds=inds;
    out.nfeats = length(selectedfeats);
    svname = fullfile('selectedfeats',sprintf('chunk%i',k));
    save(svname,'-struct','out')
    k
end

    
%%
xne = xargon('global_anls2_clust.m');
xne.argon_jobs_path = '~/xargon_submit/jobs';
% xne = xargon('global_anls2_clust_PHASERANDOMIZED.m');
xne.queue = 'NEUROSURGERY,UI,CCOM,all.q';
xne.profile = 'mid_mem';
xne.minmem = 220;

% chunk = 10;
% nchunk = 5;
% fk = 1:chunk*nchunk:length(selectedfeats);

xne.jobindices = fk;
xne.rerun_if_aborted='yes'
 xne.datafiles = 'dummy.txt';
%%
nattempts =1;
tic
while ~isempty(xne.jobindices) && nattempts<=50
    
   
    
    switch xne.status
        
        case {'running','submitted','open','waiting'}
            pause(5)
        case 'error'
            error('Status is error')
            
        case {'finished','none'}

            d = dir(fullfile(xne.subpaths.output.local,'out*.mat'));
            
            re = regexp({d.name},'out(\d*).mat','tokens','once');
            
            doneind = cellfun(@str2double,re);
            xne.jobindices = fk(~ismember(fk,doneind));
            if isempty(xne.jobindices)
                fprintf('\nFinished at %s',datestr(now))
                xne.finish()
            else
                xne.make_bash_script();
                xne.submit(true);
                fprintf('\nAttempt %i started at %0.1f s, %i remaining',nattempts,toc,length(xne.jobindices))
                nattempts = nattempts+1;
            end
        otherwise
            
            fprintf('\nStatus is %s at %s',xne.status,datestr(now))
            
    end
end


%%
% d = dir(fullfile(xne.local_save_dir,'out*.mat'));
d = dir(fullfile(xne.subpaths.output.local,'out*.mat'));

% % K3 = 0;
%  K4 = 0;

    K3 = zeros(length(selectedfeats));
    K4 = zeros(length(selectedfeats));

%%  
    badfiles ={};
    badind = [];
    
   
    
for k = 1:length(d)
% for k = redoind
    try
        ld = load(fullfile(d(k).folder,d(k).name));
        K3(ld.rows,:) = ld.K3(ld.rows,:);
        K4(ld.rows,:) = ld.K4(ld.rows,:);
        k
    catch err
        fprintf('\n%s',err.message);
        badfiles{1,end+1}=d(k).name;
        badfiles{2,end}=err.message;
        badind(end+1) = k;
    end
end
redoind = badind;

%% Randomized SVD estimation
KK = K4*K4'; % Multiply the matrices to better separate the eigenvectors
KKKK=KK*KK;  % (each multiplication squares the eigenvaue but does not change eigenvectors)
[u4,l4] = randsvd(KKKK,1000);
unorm = u4./sqrt(sum(u4.^2,2));
% [u,l] = svd(K3);

%% Write variables to file for latex

fid = fopen('techvars.tex','w');
wvar = @(name,val) fprintf(fid,'\n\def\%s{ %s }',name);

sids = unique({selectedfeats.subject});
wvar('subjectN',length(sids));
wvar('skewThresh',sprintf('0.2%f',thr));

