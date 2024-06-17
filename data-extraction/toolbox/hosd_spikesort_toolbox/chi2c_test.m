
nclust = 5;
sd = 1;
msd = 2;
dim = 3;
wishart_dof = dim + 5; %Degrees of freedom for the distribution of cluster covariances (lower -> greater variability)
max_clusts = 15;
Ntot = 5000;

seed = reseed;

rr = [0,20./Ntot+rand(1,nclust)];
Nperclust = diff(ceil(cumsum(rr)./sum(rr)*Ntot));
% Nperclust = [2000 1000 500 250 125];
% Nperclust = round(rand(1,nclust)*980)+20;
ppweight = [0 .25 .5 .75 1 nan];
isocut_threshold = [.1 .25 .5 1 2 4];
evalc_criteria = {'CalinskiHarabasz','DaviesBouldin','Silhouette'};
clusterfun =@(x,k)mogclust(x,k,'Replicates',10,'Options',struct('MaxIter',500));
% dist = 'gauss';
% dist = 'box';
dists={'box','gauss','t64','t32','t16','t8','t4','t2','t1'};

if exist('jobindex','var')
    outfile = fullfile(outputdir,sprintf('out%i.mat',jobindex));
    if exist(outfile,'file')
        return
    end
    touchdir = fullfile('~','tmp',mfilename);
    if ~exist(touchdir,'dir')
        mkdir(touchdir)
    end
    donefile = fullfile(touchdir,sprintf('out%i.done',jobindex));
    if exist(donefile,'file')
        return
    end
    touchfile = fullfile(touchdir,sprintf('out%i.touch',jobindex));
    if exist(touchfile,'file')
       d = dir(touchfile);
       if (now-d.datenum)*24<3
           return
       end
    end
    fclose(fopen(touchfile,'w'));
    
%     dists = dists(ceil(jobindex/600));
    dists = dists(mod(jobindex-1,length(dists))+1);
end

for disti = 1:length(dists)
    dist = dists{disti};
    fprintf('\nDistribution: %s',dist)
    r = [];
    trucl = [];
    for k = 1:nclust

        S = cov(randn(wishart_dof, dim ))*sd^2;
        m = randn(1,dim)*msd;

        switch dist
            case 'gauss'
                z = randn(Nperclust(k),dim);
            case 'box'
                z = zscore(rand(Nperclust(k),dim));
            otherwise
                 z = randn(Nperclust(k),dim);
                 z = z./sqrt(sum(z.^2,2));
                 z = z.*sqrt(sum(random(regexp(dist,'[^\d.]*','match','once'),str2double(regexp(dist,'[\d.]*','match','once')),Nperclust(k),dim).^2,2));
        end
        trucl = cat(1,trucl,ones(size(z,1),1)*k);

        r = cat(1,r,z*S^.5+ m);

    end

    Dmat = @(x)2*(x==x')-1-eye(length(x));
    Dtru = Dmat( trucl);
    m2lin = @(x)x(:);
    get_pwcorr = @(x)(m2lin(Dmat(x))'*Dtru(:))./sqrt(sum(abs(m2lin(Dmat(x)))).*sum(abs(Dtru(:))));
    p = @(x)x./sum(x,1);
    entr = @(x)sum(-p(x).*log2(p(x)+eps),1);
    MI = @(confus)entr(sum(confus,2))- entr(confus)*(sum(confus,1)./sum(confus(:)))';
    clear chi2out aic bic iso chz
    % for k = 1:length(ppweight)+1

    %     if k <=length(ppweight)
            [out,clval,CL]  = chi2_criterion(r,max_clusts,clusterfun,true,ppweight,0);
            gmAIC = arrayfun(@(x)x.gmAIC,clval);       
            gmBIC = arrayfun(@(x)x.gmBIC,clval);       
             [~,aic.OptimalK] = min(gmAIC);
             aic.OptimalY = clval(aic.OptimalK).gmdist.cluster(r);
             [aic.confus,aic.confusstat] = crosstab(trucl,aic.OptimalY);
             aic.MI = MI(aic.confus);
             aic.MIscore = MI(aic.confus)- MI(diag(Nperclust));
             aic.PWMcorr = get_pwcorr(aic.OptimalY);
             aic.CriterionName = 'GM_{AIC}';
             aic.CriterionValues = gmAIC;

    %           
    %         gmBIC = arrayfun(@(x)x.gmdist.BIC,clval);
            [~,bic.OptimalK] = min(gmBIC);
            bic.OptimalY = clval(bic.OptimalK).gmdist.cluster(r);
            [bic.confus,bic.confusstat] = crosstab(trucl,bic.OptimalY);
            bic.MI = MI(bic.confus);
            bic.MIscore = MI(bic.confus)- MI(diag(Nperclust));
            bic.PWMcorr = get_pwcorr(bic.OptimalY);
            bic.CriterionName = 'GM_{BIC}';
            bic.CriterionValues = gmBIC;
    %     else
    %         [out,clval]  = chi2_criterion(r,10,clusterfun,false);
    %         out.ppweight = nan;
    %         out.CriterionName = sprintf('Chi-sq. sample stats.'); 
    %     end
        for k = 1:length(ppweight)
    %         chi2out(k).ppweight = ppweight(k);
            if ~isnan(ppweight(k))
                chi2out(k).CriterionName = sprintf('Chi-sq. %0.2f',ppweight(k));
            else
                chi2out(k).CriterionName = sprintf('Chi-sq. sample');
            end
                [chi2out(k).confus,chi2out(k).confusstat] = crosstab(trucl,out.OptimalY(:,k));
            chi2out(k).PWMcorr = get_pwcorr(out.OptimalY(:,k));
            chi2out(k).MI =MI(chi2out(k).confus);
            chi2out(k).MIscore = MI(chi2out(k).confus) - MI(diag(Nperclust)); %#ok<*SAGROW>
            chi2out(k).OptimalY = out.OptimalY(:,k);
            chi2out(k).OptimalK = out.OptimalK(k);
            chi2out(k).CriterionValues = out.KSstat(:,k);
        end
    %         chi2out = out;

    % end


    for k = 1:length(evalc_criteria)
    %     eva = evalclusters(r,clusterfun,evalc_criteria{k},'Klist',1:10);
        eva = evalclusters(r,CL(:,~all(isnan(CL))),evalc_criteria{k});
        if ~isnan(eva.OptimalK)
            chz(k).OptimalY = CL(:,eva.OptimalK);
        else
            chz(k).OptimalY = nan(size(r,1),1);
        end
        chz(k).CriterionName = eva.CriterionName;
        chz(k).CriterionValues = eva.CriterionValues;

        [chz(k).confus,chz(k).confusstat] =   crosstab(trucl,chz(k).OptimalY);
        chz(k).MI = MI(chz(k).confus)';
        chz(k).MIscore = MI(chz(k).confus)' - MI(diag(Nperclust));
        chz(k).OptimalK = eva.OptimalK;
        chz(k).PWMcorr = get_pwcorr(chz(k).OptimalY);
    %     chz(k).ppweight = nan;
    end

    for k = 1:length(isocut_threshold)
        [iso(k).OptimalY,isoinfo] = isosplit5(r',struct('min_cluster_size',5,'verbose',false,'verbose_isocut',false,'isocut_threshold',isocut_threshold(k)));
        iso(k).OptimalK = length(unique(iso(k).OptimalY));
        [iso(k).confus,iso(k).confusstat] = crosstab(trucl,iso(k).OptimalY);
        iso(k).MI = MI(iso(k).confus);
        iso(k).MIscore = MI(iso(k).confus)- MI(diag(Nperclust));
        iso(k).PWMcorr = get_pwcorr(iso(k).OptimalY);
        iso(k).CriterionValues = nan;
        iso(k).CriterionName = sprintf('Isosplit5 %0.2f',isocut_threshold(k));
    end
    % iso.ppweight = nan;
    distr(disti).distribution = dist;
    distr(disti).crtests = [chi2out,chz,iso,aic,bic];
    distr(disti).r = r;
    distr(disti).trucl = trucl;
end
if exist('jobindex','var')
    save(outfile,'distr','seed','wishart_dof','Nperclust');%'chi2out','chz','seed','gmAIC','gmBIC','r','iso','aic','bic');
    fclose(fopen(donefile,'w'));
    delete(touchfile)
end

