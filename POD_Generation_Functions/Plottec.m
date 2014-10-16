function handle = Plottec(data)

% Again may need to include zones
% Also may not the data structure at all
data.x=data.xg;
data.y=data.yg;

% Again need to figure out what rmghost does
% blah = 0;
% if (size(data.x{1})) == size(data.pod),
%     blah = 1;
% end

% Determine bounds of plot
minx = min(min(data.x));
maxx = max(max(data.x));
miny = min(min(data.y));
maxy = max(max(data.y));

% TODO check what rmghost does, currently no trace of this
% if blah ==0;
%    data.y = rmghost(data.y);
%    data.x = rmghost(data.x);
% end

hold on
axis([minx maxx miny maxy])  
h = pcolor(data.x,data.y,data.pod);
colormap(jet)
shading interp  
axis equal

if min(min(data.yg)) < 0 && min(min(data.xg)) < 0 && max(max(data.xg)) > 4
    rectangle('position',[minx miny abs(minx) abs(miny)],'facecolor','k')
    rectangle('position',[4 miny (maxx-4) abs(miny)],'facecolor','k')
end
drawnow
hold off

if nargout == 1
    handle = h;
end

