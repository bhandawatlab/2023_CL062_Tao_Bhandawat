function h = filledCircle(ax,xyAll,r,matchZ,stim,varargin)
ncircles = size(xyAll,1);
if isempty(varargin)
    c = repmat([0 0 0],ncircles,1);
else
    c = varargin{1};
end
hold(ax,'on')
axis(ax,'manual')
for i = 1:ncircles
    x = xyAll(i,1);
    y = xyAll(i,2);
    th = 0:pi/50:2*pi;
    xunit = r(i) * cos(th) + x;
    yunit = r(i) * sin(th) + y;
    if matchZ(i)
        h(i) = fill(ax,xunit, yunit,c(i,:),'FaceAlpha',stim(i));
%     else
%         h(i) = plot(ax,xunit, yunit,['r:'],'LineWidth',2);
    end
%     text(ax,x+r(i),y-r(i),num2str(i),'FontSize',14,'Color','red')
end
hold(ax,'off')