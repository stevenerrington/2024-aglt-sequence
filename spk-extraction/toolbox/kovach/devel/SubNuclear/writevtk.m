
function writevtk(vtk,fn)


% Writes mesh data to an ASCII vtk file. 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/writevtk.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

%C. Kovach 2013

fid = fopen(fn,'w');

vtk.header{3} = 'ASCII';
vtk.header{2} = 'Written with writevtk.m';

for i = 1:4,
    fprintf(fid,vtk.header{i,1});
    fprintf(fid,'\n');
end

% isbinary = isequal(lower(deblank(vtk.header{3})),'binary');


fprintf(fid,'POINTS %i %s\n',vtk.npoints,vtk.fmt);
fprintf(fid,'%0.3f %0.3f %0.3f\n',vtk.vert');

fprintf(fid,'POLYGONS %i %i\n',vtk.n);
ndim = size(vtk.tri,2)-1;
str = sprintf('%i %s\n', ndim,repmat('%i ',1,ndim));

fprintf(fid,str,vtk.tri(end-ndim+1:end)'-1);

fclose(fid);

    