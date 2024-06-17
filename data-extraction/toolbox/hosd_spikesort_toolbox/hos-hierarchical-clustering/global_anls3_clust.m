% nfeats = 140318;


d = dir('~/lss/bispectral_analysis/selectedfeats/chunk*.mat');

ch1 = mod(jobindex-1,length(d))+1;
ch2 = ceil(jobindex/length(d));

out.file1 = fullfile(d(ch1).folder,d(ch1).name);
out.file2 = fullfile(d(ch2).folder,d(ch2).name);

ld1 = load(out.file1);
ld2 = load(out.file2);

for fld ={'chunk','chunkphrand'}
    fldn=fld{1};
    savename = fullfile(outputdir,sprintf('out%06i%s.mat',jobindex,fldn));
    if exist(savename,'file')
        continue
    end
    % load ~/lss/bispectral_analysis/selectedfeats.mat



    % af1 = allfeats;
    % af1 = selectedfeats;
    % 
    % af1 = allfeats([allfeats.component]==1);
    % 
    % mx = [af1.featfiltmax];
    % [srt,srti] = sort(mx,'descend');
    % selectedfeats = af1(srti);

    % compno = [selectedfeats.component];
    F1 = [ld1.(fldn).feature];
    F1(isnan(F1)) = 0;
    F2 =[ld2.(fldn).feature];
    F2(isnan(F2)) = 0;
    PDF1 = [ld1.(fldn).pdfilt];
    PDF1(isnan(PDF1)) = 0;
    PDF2 = [ld2.(fldn).pdfilt];
    PDF2(isnan(PDF2)) = 0;

    F1 = ifft(fft(F1).*abs(fft(PDF1)));
    F1 = zscore(F1);
    F2 = ifft(fft(F2).*abs(fft(PDF2)));
    F2 = zscore(F2);

    ftF1 = fft(F1);
    ftF2 = fft(F2);
    % ftPDF = fft(PDF);
    FF1 = ifft(abs(ftF1).^2);
    FF2 = ifft(abs(ftF2).^2);
    % FF = ifft(ftF.*ftPDF);
    s13 = sum(FF1.^3);
    s14 = sum(FF1.^4);
    mxFF1 = max(FF1);
    s23 = sum(FF2.^3);
    s24 = sum(FF2.^4);
    mxFF2 = max(FF2);

    if exist('jobindex','var')
        out.K3 = sparse(ld1.nfeats,ld2.nfeats,size(F1,2).*size(F2,2));
        out.K4 = out.K3;
        out.M3 = out.K3;
        out.M4 = out.K3;
        out.Kpeakcc = out.K3;
    else

        K3 = zeros(size(F,2));
        K4 = zeros(size(F,2));
        Kpeakcc = K3;
    end
    % ftPDFperm = permute(ftPDF,[1 3 2]);
    % ftPDFperm = permute(conj(ftF),[1 3 2]);

    tic
%     rows = [];
%     for chnki = 1:nchunks

    k3 = cumdist(F1,F2,3);
    m3 = momdist(F1,F2,3);
    out.K3(ld1.inds,ld2.inds) = real(k3);
    out.M3(ld1.inds,ld2.inds) = real(m3);
    [k4,dsgn,pkcc] = cumdist(F1,F2,4);
    m4 = momdist(F1,F2,4);
    out.K4(ld1.inds,ld2.inds) =real(k4);
    out.M4(ld1.inds,ld2.inds) = real(m4);
    out.Kpeakcc(ld1.inds,ld2.inds) = pkcc;

    out.rows = ld1.inds; 
    out.cols = ld2.inds; 

%     k3(isnan(k3)) = 0;
%     k4(isnan(k4)) = 0;
%     m3(isnan(k3)) = 0;
%     m4(isnan(k4)) = 0;
%     [u,l,v] = svd(k3);
%     [u,l,v] = svd(k4);
     
    
    toc
     if exist(savename,'file')
        return
     end
    try
        fid = fopen([mfilename,'.m'],'r');
        out.COM = fread(fid,'uchar=>char')';
    catch
        out.COM = '';
    end


    save(savename,'-struct','out');

end