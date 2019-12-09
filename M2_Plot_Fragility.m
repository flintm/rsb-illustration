function M2_Plot_Fragility(x_lim,Sa,P_DS,lty,color,DS,name,struc, ylabT, xlabT, saveFig)
Lwds = [0.5, 1, 1.5]; % PP of 1, 1.5, 2

gcf = figure('Color',[1 1 1]);
set(gcf,'units','inches','position',[1 1 1.5 1.5],'PaperUnits', 'Inches');
gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
            'FontSize',10,...
    'FontName','Arial',...
    'Linewidth', 1,...
    'TickLength', [0.02 0.035],...
    'XLim',x_lim);%,... % 'YLim',y_lim,
   % 'XTick',[x_lim(1)+0.5,mean(x_lim),x_lim(2)-0.5]);%,...
   % 'YTick',[y_lim(1)+0.5,mean(y_lim),y_lim(end)-0.5]);
set(gca, 'units','inches','position',[0.45 0.45 0.9 0.9])
box(gca,'on');
hold(gca,'all');

if ~xlabT
    set(gca,'XTickLabel',[]);
end
if ~ylabT
 set(gca,'YTickLabel',[]);
end

for i = 1:3
    plot(Sa,P_DS(:,i),'LineStyle',lty,'Color',color,'Linewidth',Lwds(i));
end

if xlabT
    xlabel('Sa [g]','FontSize',10,'FontName','Arial')
end
if ylabT
    ylabel('P(ds > DS|IM)','FontSize',10,'FontName','Arial')
end

set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')

text(x_lim(2)*0.95,0.1,DS,'HorizontalAlignment','right')

print('-depsc', [pwd '/Figs/eps/frag_',DS,'_',name,'_',struc,'.eps']);
if saveFig==true
    savefig([pwd '/Figs/fig/frag_',DS,'_',name,'_',struc,'.fig'])
end
end

