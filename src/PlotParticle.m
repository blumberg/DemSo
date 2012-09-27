ang = linspace(-pi,pi,100);
x = cos(ang)*Rx + Center(1);
y = sin(ang)*Ry + Center(2);
plot(x,y,'-k')
hold on
plot(pos(1),pos(2),'o','LineWidth',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerSize',60);
set(gcf,'Renderer', 'OpenGL');
set(gcf,'Toolbar', 'none');
xlim([0,10]);
ylim([0,10]);
grid off
box off
axis off
hold off
set(gca,'Position',[0.0,0.0,1,1])
set(gcf,'Renderer', 'OpenGL');
set(gcf,'Toolbar', 'none');

drawnow
