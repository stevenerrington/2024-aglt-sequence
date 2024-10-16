
function vtk = readvtk(fn)

%
% Reads vtk files
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/readvtk.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2013

fid = fopen(fn,'r');

for i = 1:4,
    vtk.header{i,1} = fgetl(fid);
end

isbinary = isequal(lower(deblank(vtk.header{3})),'binary');


currfield = 'xxx';
doread = true;
while ~feof(fid) 
    
    l = fgetl(fid);
    
    if feof(fid) || isequal(l,-1)
        break
    end
        
    if ~isbinary

        scl = sscanf(l,'%f');

        if isempty(scl)
            scl = sscanf(l,'%s',1);
            switch scl
                case 'POINTS'
                    currfield = 'vert';
                    vtk.npoints = sscanf(l,'%*s %i %*s');
                    vtk.fmt  = char(sscanf(l,'%*s %*i %s')');
                    fun=@(x)x;
                case 'POLYGONS'
                    vtk.n = sscanf(l,'%*s %i %i')';
                    currfield = 'tri';
                    fun = @(x)x(2:end)+1;
                otherwise
                    currfield = scl;
            end
                vtk.(currfield) = [];
        else
            vtk.(currfield)(end+1,:) =fun(scl);
        end
    else
            scl = sscanf(l,'%s',1);
    
            switch scl
                case 'POINTS'
                    currfield = 'vert';
                    vtk.npoints = sscanf(l,'%*s %i %*s');
                    fmt  = char(sscanf(l,'%*s %*i %s')');
                    fun = @(x)reshape(x,3,length(x)/3)';
                    readn = vtk.npoints*3;
                case 'POLYGONS'
                    n = sscanf(l,'%*s %i %i')';
                    vtk.ntess = n(1);
                    currfield = 'tri';
                    fmt = 'int32';
                    readn = n(2);
                    rmcol = @(x)x(:,2:end);
                    fun = @(x) rmcol(reshape(x,4,length(x)/4)')+1;
                case 'FIELD'
%                     currfield = 'field';
                    fldn = sscanf(l,'%*s %*s %i')';
                    fldi = 1;
                    currfield = char(sscanf(l,'%*s %s %*i')');
                    readn = 0;
                    fmt = 'int';
                    fun=@(x)[];
                otherwise
                    if fldi <=fldn %#ok<BDSCI>
                        n2 = sscanf(l,'%*s  %i %i %*s')';
                        fnm = sscanf(l,'%s',1); 
                        fmt = char(sscanf(l,'%*s  %*i %*i %s')');
                        readn = prod(n2);
                        fun = @(x) cat(2,vtk.(currfield),struct('label',fnm,'data',reshape(x,n2)));
                        fldi = fldi+1;
                    else
                        error('unrecognized data type')
                    end
            end
            
             vtk.(currfield) = fun(fread(fid,readn,fmt));
    end

end
fclose(fid)
    