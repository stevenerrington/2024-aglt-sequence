 function plotupdate(me,axs,updatesis)

 %%% Update Plots
 
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/@volumeview/plotupdate.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

 me.cleanup();
 if nargin < 3
     updatesis = true;
 end

 if nargin < 2 || isempty(axs)
    axs = 1:length(me.plotax); 
 end
 crp = me.current_point;

for k = 1:length(me.modules)
    if isa(me.modobjects.(me.modules{k}),'SubModules.module')
      me.modobjects.(me.modules{k}).updateview(); 
    end
end
 for k = 1:length(me.volumes)
     if ~me.volumes(k).show
         set(me.volumes(k).ploth,'visible','off')
         continue
     else
        for i = axs(:)'
             nv =me.plotax(i).normvec;
             rot = @(x)rot90(x,me.plotax(i).rot);
             if get(me.plotax(i).handles.transpose,'value')
                trfunm = @(x)rot(x);                                    
             else
                 trfunm = @(x)rot(x');               
             end
             if i<=3 %First thre axes are the standard orthogonal image planes
               dim = (1:3)*nv';


                 xix = {':',':',':'};
                 xix{dim} = round(crp(dim)+1);

                 x = squeeze(me.Vol(xix{:}));


             else


                    vc = size(me.current.volumes.image.Data)/2;
                   p0 = ((crp-vc)*nv')*nv+vc;
                  [x,~,trmat] = me.makeslice(nv,p0);
                  me.plotax(i).basetrmat=trmat;
                  
%                    me.rtrmat(me.plotax(i));
             end
             x = trfunm((x-me.volumes(k).intensity_range(1))/diff(me.volumes(k).intensity_range));
             x(x>1) = 1;
             x(x<0) = 0;

            me.volumes(k).ploth = me.volumes(k).ploth(ishandle(me.volumes(k).ploth));
                set(me.volumes(k).ploth(i),'cdata',x(:,:,[1 1 1]),'visible','on','xdata',[1 size(x,2)],'ydata',[1 size(x,1)]);                          
            
            
    if get(me.fixSisterAx,'value') && updatesis
             sisax = me.plotax(i).sisters;
         if ~isempty(sisax)
             sisax(~isvalid(sisax))=[];
             sisax(~ishandle([sisax.h])) = [];
             for kk = 1:length(sisax)
                    axis(sisax(kk).h,axis(me.plotax(i).h));

             end
             me.plotax(i).sisters=sisax;
         end
    end
    
     crk = me.plotax(i).handles.crosshairs;   
     crk = crk(ishandle(crk));
     if get(me.showpt,'value')
         othax = find([me.plotax.h]~=me.plotax(i).h);
         for kk = 1:length(othax)
             ko = othax(kk);
              axl=axis(me.plotax(ko).h);
              cpax = me.plotax(ko).vol2ax(me.current_point)';
              xc = [axl(1:2)',cpax([ 2 2]);cpax([1 1]),axl(3:4)'];
            nv2 = me.plotax(ko).normvec;
            vx = me.plotax(ko).ax2vol(xc);
            q=cross(nv,nv2);
            q=q./norm(q);
            qcx= ((vx-repmat(crp,size(vx,1),1))*q');
            [~,mni] = min(qcx);
            [~,mxi] = max(qcx);
            vxpr = qcx([mni mxi],:)*q+repmat(crp,2,1); 
             crlim = me.plotax(i).vol2ax(vxpr);
             set(crk(kk),'xdata',crlim(:,1),'ydata',crlim(:,2),'visible','on','zdata',[.5 .5])
         end
     else
         set(crk(:),'visible','off')
     end 

        end
    end


 end
 me.plotmeshes;
 arrayfun(@(a,b) set(a,'String',round(b)),me.pointbox,round(crp));
 arrayfun(@(a,b) set(a,'String',b),me.coordbox,round(me.current.transforms(1).tr(crp)));
 me.annotationUpdate();
       plhs = cat(1,me.meshes.ploth,me.points.ploth);
       shown = cat(1,me.meshes.show,me.points.show)==1;
     if ~isempty(plhs(shown))  
         lbls = {me.meshes.label,me.points.label};
         if get(me.showlegendh,'value') ==1
            l =  legend(plhs(shown,3),lbls(shown),'location','northOutside');
         elseif ~isempty(me.legendh) && ishandle(me.legendh)
             delete(me.legendh);
         end
             
             
             
        pos = get(l,'position');
        me.legendh = l;
        set(l,'units','normalized','location','none','position',pos + [0 .1 0 0]) 
     else
         legend off
     end
     