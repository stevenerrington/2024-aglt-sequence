function out = readurl_binary(url)

% Retrieve binary data from a a url
%
u = java.net.URL(url);
uc = u.openConnection();
out.contentType = char(uc.getContentType());
out.contentLength = uc.getContentLength();
raw = uc.getInputStream();
in = java.io.BufferedInputStream(raw);

out.data = zeros(1,out.contentLength,'uint8');
% tic
k = 1;
dk = 1;
dispn = round((1:100)/100*out.contentLength);
h = waitbar(0,'Retrieving binary data...slowly.');
n = in.read;
while n~=-1
    out.data(k) = n;
    n = in.read;
   
    if k==dispn(dk)
        waitbar(k/out.contentLength,h)
        dk=dk+1;
    end
     k = k+1;
%     fwrite(fid,in.read);
end
delete(h)
% toc