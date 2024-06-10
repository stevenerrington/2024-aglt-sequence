function templ = make_template(urldata)


% tmpl = make_template(urldata)
%
% This function takes data from the protocol url and generates a template
% structure.
%
%
% tmpl contains template data for each channel group.
%
% see also parse_block_url

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/make_template.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

writefile = true;
if nargin < 1
    urldata = parse_block_url;
elseif ischar(urldata)
    urldata = parse_block_url(urldata);
end

lozch = urldata.lozchannels;
labels = {lozch.label};
labels(cellfun(@(x)~ischar(x),labels))={'ZZZ_EMPTY'};

[unqlbl,~,lbli]=unique(labels);


z = 0;
dput = [0 0];

if writefile
    fn = sprintf('template%s.m',urldata.subject_id);
    [fn,pth] = uiputfile(fn,'Please choose a template Name');
    if ~isempty(unqlbl)
        fid = fopen(fullfile(pth,fn),'w');
    else
        error('No electrodes identified.')
    end
    fprintf(fid,'function tmp=%s\n\n',strtok(fn,'.'));
    fprintf(fid,'\n\n%% This is a template file generated for %s.\n',urldata.subject_id);
    fprintf(fid,'%%   tmp = tmp=%s\n\n%%  tmp is a structure containing a template of electrode\n%% positions and interelectrode distances.',strtok(fn,'.'));
    fprintf(fid,'\n\n\nindx=0;');
    fprintf(fid,'\nprurl=''%s'';',urldata.url);
    fprintf(fid,'\n\nz=%g;',z);
    % fprintf(fid,'\ntmp(indx).subject=''%s'';',templ.subject);
    % fprintf(fid,'\ntmp(indx).protocol_url=prurl;');

end

 templ.subject = urldata.subject_id;
templ.protocol_url= urldata.url;
   

for i = 1:length(unqlbl)
    
    if strcmp(unqlbl{i},'ZZZ_EMPTY')
        continue;
    end
    geti = lbli==i;
    ncontact = sum(geti);
    
    group(i).label = unqlbl{i};
    group(i).channels= lozch(geti);
    group(i).ncontact = ncontact;
    mktmp = true;
    switch ncontact
        case 4 %%% 4 contact strip
            dx = 10; %%% inter electrode distance in mm
            dimensions = [1 4];
        case 10 %%% parahipocampal
            dx = 10;
            dimensions = [1 10];
        case 6  %%% some depth 
            dx = 10;
           dimensions = [1 6];

        case 96  %%% 96 contact grid
            dx = 5;
           dimensions = [8 12];
         
        case 32  %%% 32 contact frontal grid
            dx = 10;
           dimensions = [4 8];
            
        case 16  %%% 16 contact  grid
            dx = 5;
            dimensions = [2 8];
            
        case {23 27} %%% Frontal strips
            
            dx = 10;
            [Y,X] = meshgrid((1:6),(1:6)); %temporal pole will be modeled as separate strips
            keep = (X<2|Y<6) & (X<3|Y<5); 
            Y = Y+dput(2)/dx;
            X = X+dput(1)/dx;
            XY = [X(keep),Y(keep)]*dx;
            group(i).Ds = dfun(XY,XY);
            group(i).XY(:,[1 3]) = XY;
            group(i).XY(:,2) = z;
            group(i).dimensions = size(keep);
            group(i).dx = dx;
            group(i).dput = dput;
            mktmp = false;
            if writefile
                write_file(fid,group(i),true);
            end

        otherwise
            if mod(ncontact,8)==0
                dx =5; %% other grids/strips assume 5mm
                dimensions = [ncontact/8,8];
            else
                warndlg(sprintf('Unexpected number of contacts, %i, for %s. Please check.',ncontact,unqlbl{i}))
            end
    end
    
    if mktmp
        group(i).dx = dx; %#ok<*AGROW>
        group(i).dimensions = dimensions;
        [Y,X] = meshgrid(dput(1)/dx+(1:group(i).dimensions(2)),(1:group(i).dimensions(1))+dput(2)/dx); %16-contact grid
        XY = [X(:),Y(:)]*dx;
        group(i).Ds = dfun(XY,XY);
        group(i).XY(:,[1 3]) = XY;
        group(i).XY(:,2) = z;
        group(i).dput = dput;
        if writefile
            write_file(fid,group(i),false);
        end
%         hold on
%         plot(group(i).XY(:,1),group(i).XY(:,3),'.'); pause
    end
    if mod(i,8)==0;
        mx = max(cat(1,group.XY));
        dput(1) = mx(3)+ 20;
        dput(2)=0;
    else
        dput(2)=dput(2)+4*dx*dimensions(1);
    end
            
end

n = [group.ncontact];
[srt,srti] = sort(n,'descend');
group = group(srti);

templ.group=group;
T = blkdiag(group.Ds);
T(T==0) = -1;
T = T.*(1-eye(length(T)));
T(end+1,:) = -1;
T(:,end+1) = -1;
templ.T=T;

if writefile
    fprintf(fid,'\nn = [group.ncontact];');
    fprintf(fid,'\n[srt,srti] = sort(n,''descend'');');
    fprintf(fid,'\ngroup = group(srti);');
    fprintf(fid,'\ntmp.group=group;');

    fprintf(fid,'\n\n\n%% Template with all channels');
    fprintf(fid,'\nT = blkdiag(group.Ds);');
    fprintf(fid,'\nT(T==0) = -1;');
    fprintf(fid,'\nT = T.*(1-eye(length(T)));');
    fprintf(fid,'\nT(end+1,:) = -1;');
    fprintf(fid,'\nT(:,end+1) = -1;');
    fprintf(fid,'\ntmp.T=T;');

    fclose(fid);
end
%%%%%%

function write_file(fid,templ,tpole)

if nargin < 3
    tpole= false;
end
    
% Write the file
fprintf(fid,'\n\n\n%%%%%%%%%% %s %%%%%%%%%%\n',templ.label);
fprintf(fid,'\ndput=[%g %g];',templ.dput);
fprintf(fid,'\nindx=indx+1;');
fprintf(fid,'\ngroup(indx).ncontact=%i;',templ.ncontact);
fprintf(fid,'\ngroup(indx).label=''%s'';',templ.label);
% fprintf(fid,'\ngroup(indx).subject=''%s'';',templ.subject);
% fprintf(fid,'\ngroup(indx).protocol_url=prurl;');
fprintf(fid,'\ngroup(indx).channels=[];');
fprintf(fid,'\ndx=%g;',templ.dx);
fprintf(fid,'\ngroup(indx).dx=dx;');
fprintf(fid,'\ngroup(indx).dimensions=[%i %i];',templ.dimensions);
if tpole
    fprintf(fid,'\n[Y,X] = meshgrid((1:6),(1:6)); %%temporal pole will be modeled as separate strips');
    fprintf(fid,'\nkeep = (X<2|Y<6) & (X<3|Y<5);'); 
    fprintf(fid,'\nY = Y+dput(2)/dx;');
    fprintf(fid,'\nX = X+dput(1)/dx;');
    fprintf(fid,'XY = [X(keep),Y(keep)]*dx');
        
else
    fprintf(fid,'\n[Y,X] = meshgrid(%i/dx+(1:group(indx).dimensions(2)),(1:group(indx).dimensions(1))+%i/dx);',templ.dput);
    fprintf(fid,'\nXY = [X(:),Y(:)]*group(indx).dx;');
end
fprintf(fid,'\ngroup(indx).Ds = dfun(XY,XY);');
fprintf(fid,'\ngroup(indx).XY(:,[1 3]) = XY;');
fprintf(fid,'\ngroup(indx).XY(:,2) = z;');
fprintf(fid,'\ngroup(indx).dput = dput;');

% fprintf(fid,'\n\nif mod(i,4)==0;');
% fprintf(fid,'\n\tmx = max(cat(1,group(indx).dimensions));');
% fprintf(fid,'\n\tdput(2) = mx(2)*5 +10;');
% fprintf(fid,'\n\tdput(1)=0;');
% fprintf(fid,'\nelse');
% fprintf(fid,'\n\tdput(1)=dput(1)+dx*dimensions(1);');
% fprintf(fid,'\nend\n\n\n');
% 


