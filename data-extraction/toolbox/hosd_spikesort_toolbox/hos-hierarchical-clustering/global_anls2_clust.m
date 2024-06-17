
savename = fullfile(outputdir,sprintf('out%06i.mat',jobindex));
if exist(savename,'file')
    return
end
% load ~/lss/bispectral_analysis/selectedfeats.mat
load ~/bispectral_analysis2/selectedfeats.mat

chunk = 10;
nchunks = 5;

% af1 = allfeats;
% af1 = selectedfeats;
% 
% af1 = allfeats([allfeats.component]==1);
% 
% mx = [af1.featfiltmax];
% [srt,srti] = sort(mx,'descend');
% selectedfeats = af1(srti);

compno = [selectedfeats.component];
F = [selectedfeats.feature];
F(isnan(F)) = 0;
PDF = [selectedfeats.pdfilt];
PDF(isnan(PDF)) = 0;

F = ifft(fft(F).*abs(fft(PDF)));
F = zscore(F);

ftF = fft(F);
% ftPDF = fft(PDF);
FF = ifft(abs(ftF).^2);
% FF = ifft(ftF.*ftPDF);
s3 = sum(FF.^3);
s4 = sum(FF.^4);
mxFF = max(FF);

if exist('jobindex','var')
    K3 = sparse(size(F,2),size(F,2),nchunks*chunk*size(F,2));
    K4 = sparse(size(F,2),size(F,2),nchunks*chunk*size(F,2));
    Kpeakcc = K3;
else
    
    K3 = zeros(size(F,2));
    K4 = zeros(size(F,2));
    Kpeakcc = K3;
end
% ftPDFperm = permute(ftPDF,[1 3 2]);
% ftPDFperm = permute(conj(ftF),[1 3 2]);

tic
rows = [];
for chnki = 1:nchunks
   
    fk = jobindex + chunk*(chnki-1);
    inds = fk:min(fk+chunk-1,size(F,2));
    ftFsub = ftF(:,inds);
%     krind = kron(ones(1,size(ftFsub,2)),1:size(ftPDF,2));
    
%     fF =   ifft(ftFsub.*ftPDFperm);
%     Q = repmat(ftFsub,1,size(ftPDF,2)).*ftPDF(:,krind);
%     fF = ifft(ftF.*conj(ftF(:,fk)));
%     fF = ifft(ftF.*ftPDF(:,fk)); % This implicitly normalizes by SNR.
    k3 = cumdist(F,F(:,inds),3);
    K3(inds,:) = k3';
    [k4,dsgn,pkcc] = cumdist(F,F(:,inds),4);
    K4(inds,:) =k4';
    Kpeakcc(inds,:) = pkcc';
   
    rows = [rows,inds]; %#ok<AGROW>
    fk %#ok<NOPTS>
    toc
     if exist(savename,'file')
        return
    end
end

save(savename,'K3','K4','Kpeakcc','rows');

