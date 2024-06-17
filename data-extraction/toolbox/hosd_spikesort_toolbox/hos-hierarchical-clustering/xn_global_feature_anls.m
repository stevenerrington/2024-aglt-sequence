
% [res,out] = system('find ~/lss/bispectral_analysis/ -maxdepth 1 -regex .*[0-9]*-[0-9]*');
[res,out] = system('find ~/bispectral_analysis2/the_big_analysis/ -maxdepth 1 -regex .*[0-9]*-[0-9]*');

re = regexp(out,'[^\n]*','match');

[res,out] = system(sprintf('find the_big_analysis/  -maxdepth 2 -regex .*hos.*[.]mat'));
nfiles = length(regexp(out,'[^\n]*'));

getatlas = {'BN246','KN_Anatomical_region','Yeo_7_Network','Yeo_17_Network'};
mw = mwapi2;
 wref = 0:.25:150;
 tref = -2:.002:2; 
%  url = 'https://saccade.neurosurgery.uiowa.edu/labwiki/index.php/Subject/';
allfeats = [];
nproc = 0;
tic

 %%
 
inputdat=[];
subdat = struct();
for k = 1:length(re)
    
    [res,out] = system(sprintf('find %s -maxdepth 3 -regex .*hos.*[.]mat',re{k}));
  
    re2 = regexp(out,'[^\n]*','match');
    re2 = regexprep(re2,'/home/kovachc','~');
%     re3 = regexp(re2,'/([^/]*)/[^/]*$','tokens','once');
%     re3  = [re3{:}];
%     [unq,unqi] = unique(re3);
    blk = regexp(re{k},'(\d\d\d)-(\d\d\d)','tokens','once');
   sid = blk{1};
   fld = sprintf('s%s',sid);
   if isfield(subdat,fld)
       cdats = subdat.(fld);
   else
       cdats = mw.askargs({'Subject',sid,'Category','Electrode Contact'},[{'Contact Number','At CIT168toMNI X','At CIT168toMNI Y','At CIT168toMNI Z','DKT label','Destrieux label'},strcat(getatlas,' label')],{'link','none'});
       subdat.(fld) = cdats;
   end
   
   if isempty(cdats)
       continue
   end
   
     inputdat = [inputdat,struct('sid',sid,'cdats',cdats,'blk',{blk},'getatlas',{getatlas},'tref',tref,'wref',wref,'re2',{re2})];
   fprintf('\n%i of %i',k,length(re))
end

save contact_data inputdat

%%
xne = xargon('global_feature_anls_clust.m');
% xne = xargon('global_feature_anls_clust_PHASERANDOMIZED.m');
xne.queue = 'NEUROSURGERY,UI-HM,CCOM,all.q';
% xne.queue = 'UI,CCOM,all.q';
xne.profile = 'high_mem';
xne.nslots = 4;
xne.nparallel = length(inputdat);
xne.datafiles = 'contact_data.mat';
xne.submit;
