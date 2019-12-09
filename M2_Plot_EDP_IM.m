function M2_Plot_EDP_IM(x_lim,y_lim,Sa,EDP,Sa_c,color,name,struc, ylabT,xlabT,saveFig)

gcf = figure('Color',[1 1 1]);
    set(gcf, 'units','inches','position',[1 1 1.77 1.5],'PaperUnits', 'Inches');
    gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
                'FontSize',10,'Color','none',...
        'FontName','Arial',...
        'Linewidth', 1,...
        'TickLength', [0.02 0.035],...
        'YLim',y_lim,'XLim',x_lim,'XColor',[0,0,0],...
        'Xscale','log','Yscale','log','YColor',[0,0,0],...
        'XTick',[.02, 0.1, 1],...
        'YTick',[0.001,0.01,0.1]);
    box(gca,'on');
    hold(gca,'all');
if ~xlabT
    set(gca,'XTickLabel',[]);
else
    lb = get(gca,'XTick');
    for i=1:size(lb,2)
       lbn{i} = num2str(lb(i));
    end
    set(gca,'XTickLabel',lbn);
end
if ~ylabT
 set(gca,'YTickLabel',[]);
else
    lb = get(gca,'YTick');
    for i=1:size(lb,2)
       lbn{i} = num2str(lb(i));
    end
    set(gca,'YTickLabel',lbn);
end

plot(gca,[x_lim(1),x_lim(2)],[(Sa_c),(Sa_c)],'--k')

if strcmp(struc,'conc')
    scatter(gca,(Sa),(EDP(:,1)),'o','MarkerEdgeColor',color,'LineWidth',0.5);
else
    scatter(gca,(Sa),(EDP(:,1)),'s','MarkerEdgeColor',color,'LineWidth',0.5);
end
scatter(gca,(Sa),(EDP(:,2)),'x','MarkerEdgeColor',color,'LineWidth',0.5);
scatter(gca,(Sa),(EDP(:,3)),'+','MarkerEdgeColor',color,'LineWidth',0.5);

if xlabT
    xlabel('Sa','FontSize',10,'FontName','Arial')
end
if ylabT
    ylabel('Drift','FontSize',10,'FontName','Arial')
end

text(x_lim(1)*2,Sa_c*1.7,'C')
text(x_lim(1)*2,Sa_c*0.6,'NC')
  
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')

print('-depsc', [pwd,'/Figs/eps/',name,'_',struc,'.eps']);
if saveFig==true
    savefig([pwd '/Figs/fig/',name,'_',struc,'.fig'])
end
end

