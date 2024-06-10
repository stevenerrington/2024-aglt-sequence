function [IM,hdat] = loadimg(fname)

fnoext = strtok(fname,'.');
hdat = readhdr([fnoext,'.hdr']);

% datafmts = {'ubit1','uchar','int16','int32','single','double','double[2]','double[3]','char','uint16','uint32','int32','uint64','float128','double[2]','float128[2]'};
datafmts = {'ubit1','uchar','int16','int32','single','double','double[2]','double[3]'};
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/loadimg.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


fid = fopen([fnoext,'.img']);

if hdat.datatype == 0
   warning('Unknown data type')
   lbl = strcat(cellfun(@num2str,num2cell(1:length(datafmts)),'uniformoutput',false),{'. '},datafmts);
   inp = input(['\nAssume which type?',sprintf('\n%s',lbl{:}),'\n']);
   format = datafmts{inp};
else
    
    format = datafmts{log(double(hdat.datatype))./log(2)+1};
end

img = fread(fid,format);
fclose(fid)

imgdim = double([hdat.dim1, hdat.dim2, hdat.dim3]);

IM = reshape(double(img),imgdim);


