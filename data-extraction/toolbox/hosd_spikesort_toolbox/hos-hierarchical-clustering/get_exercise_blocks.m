
inputdir = 'the_big_analysis';
ld = load('exercise_bocks.mat');
blks = ld.q;

for k = 1:length(blks)
% k    
    d = dir(fullfile(inputdir,blks(k).main,'*hos.mat'));
    fprintf('\n%s chan ',blks(k).main);
    nfp =0;
    for kk = 1:length(d)    
         ld = load(fullfile(d(kk).folder,d(kk).name));
         ld.sub = ld.block.subject;
         ld.blockid = ld.block.block;
         ld.proto = ld.block.subprotocol;
         ld.contact= ld.chan.contact;
         
        if kk ==1 && k ==1
            lds = ld;
        else        
            lds(end+1) = ld;
        end
        nfp = fprintf([repmat('\b',1,nfp),'%i'],kk)-nfp;
        
    end
end


hoss = cat(1,lds.hos);
%%
fsref = 500;
tref = -2:1/fsref:2;
F = [];
PDF = [];
skew = [];
kurt = [];
for k = 1:length(lds)
    F(:,k,:) = ftresamp([lds(k).hos.feature],ld.fs,fsref,length(tref));
    PDF(:,k,:) = ftresamp([lds(k).hos.filterfun],ld.fs,fsref,length(tref));
%     skew(k,:) = skewness(lds(k).hos.xfilt(ld.dat));
%     kurt(k,:) = kurtosis(lds(k).hos.xfilt(ld.dat));
    
    k
    if k==1
        F(:,length(lds),:) = 0;
        PDF(:,length(lds),:) = 0;
    end
end
Fnorm = zscore(ifft(fft(F).*abs(fft(PDF))));
Ffilt = ifft(fft(F).*fft(PDF));

tab = readtable('Exercise_data_contact_labels.xls');

ar = arrayfun(@(x)find(contains(tab.Subject,x.sub) & tab.ContactNumber==x.contact),lds); 

isekg = contains(tab.BelongsToContactGroup(ar),'ekg','IgnoreCase',true);
issz = contains(tab.KNAnatomicalRegionLabel(ar),'sz','IgnoreCase',true);
ishipp = contains(tab.KNAnatomicalRegionLabel(ar),'hipp','IgnoreCase',true);
sub = tab.Subject(ar);

[srt, srti] = sort(squeeze(max(Ffilt)),'descend');

[unq,~,unqi] = unique({lds.blockid});
skew = zeros(length(lds),3);
for k = 1:length(unq)
    X = [lds(unqi==k).dat];
    X(isnan(X)) = 0;
    ff=PDF(:,unqi==k,:);
    ff(length(X),:,:) = 0;
    Xf = real(ifft(fft(X).*fft(ff)));
    skew(unqi==k,:) = skewness(Xf);
    
    k
end
    
[srt,srti] = sort(skew,'descend');

%%
anova1(skew(:,1),ishipp + 2*issz);
tab2 = tab(ar,:);
tab2.Seizure = issz;
tab2.Hippocampus = ishipp;

% tab2.Protocol= {lds.proto}';
blkprotos = [repmat({'Pre resting'},17,1),repmat({'Exercise'},17,1),repmat({'Post resting'},17,1)]';
tab2.Protocol = blkprotos(unqi);
tab2.Contact = ar';
[subs,~,subi] = unique({lds.sub});
[prots,~,proti] = unique({lds.proto});
    

ft = fitlme(tab2,'skew~Protocol*(Seizure + Hippocampus)  + (1|Contact) + (Seizure + Hippocampus|Subject) ');
%%
figure, 
bh = bar(ft.fixedEffects)
hold on, 
plot([1 1]'*(1:length(ft.fixedEffects)),ft.coefCI','k');
set(gca,'xtick',1:length(ft.fixedEffects),'xticklabel',regexprep(ft.CoefficientNames,'_','\\_'),'xticklabelrotation',30)
