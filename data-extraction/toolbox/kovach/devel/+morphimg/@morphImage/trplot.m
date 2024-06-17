function trplot(morphim,bh,ev,hd) %#ok<INUSD>

h = get(bh,'parent');
% trh = getappdata(h,'trh');
trh = morphim.trh;
if ~isempty(trh)&&all(ishandle(trh(:))), delete(trh), end
if get(getappdata(h,'buttonh'),'Value')==0,
    return,
end
% VX = morphim.VX;
% VY = morphim.VY;
% axs = getappdata(h,'axes');
axs = morphim.axes;
% wgt = morphim.tesswgt;


% if isempty(VX), return, end
VM = morphim.images.weightedsum(morphim.tesswgt);
D=DelaunayTri(VM);

for i = 1:length(morphim.images)
    axes(axs(i)) %#ok<MAXES>
    morphim.trh(:,i) = triplot(D.Triangulation,morphim.images(i).vert(:,1),morphim.images(i).vert(:,2));
end
% setappdata(h,'trh',trh)



