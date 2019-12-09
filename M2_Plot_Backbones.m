%% Setting backbone parameter values for 9 cases
close all
struc = 'steel';
if strcmp(struc,'conc')
    Lstyle='-';
    Mstyle='o';
    
    %  critical damping for concrete
    dampRatio = 0.05;
    %  different SDOFs' backbone parameters
    %  Displacements are normalized by structure height
    %  Forces are base shear normalized by weight
    %  yield displacement
    dy_list=[0.005 0.005 0.005 0.0075 0.01 0.005 0.005 0.0075 0.01 ];
    %  capping displacement
    dc_list=[0.046 0.046 0.046 0.049 0.051 0.067 0.087 0.069 0.092];
    %  ultimate displacement
    du_list=[0.1 0.13 0.15 0.10 0.10 0.12 0.14 0.15 0.2];
    %  yeild force list
    Fy_list=[0.147 0.147 0.147 0.221 0.294 0.147 0.147 0.221 0.294];
    %  capping force list
    Fc_list=[0.133 0.133 0.133 0.2 0.26 0.133 0.133 0.2 0.27];
    %   K [expr (2*$pi/1.16)**2]
else
    Lstyle='--';
    Mstyle='s';
    dampRatio = 0.02;
    dy_list=[0.009 0.009 0.009 0.0135 0.018 0.009 0.009 0.0135 0.018];
    dc_list=[0.03 0.03 0.03 0.035 0.039 0.041 0.051 0.045 0.06];
    du_list=[0.08 0.11 0.13 0.08 0.09 0.09 0.1 0.12 0.16];
    Fy_list=[0.2 0.2 0.2 0.3 0.4 0.2 0.2 0.3 0.4];
    Fc_list=[0.19 0.19 0.19 0.29 0.38 0.19 0.19 0.29 0.38 ];
    %	 K [expr (2*$pi/1.32)**2]
end

%   this factor accounts for the different between the location of literature and our study
ScaleFactor = 0.46;
dy_list = dy_list*ScaleFactor;
dc_list = dc_list*ScaleFactor;
du_list = du_list*ScaleFactor;
Fy_list = Fy_list*ScaleFactor;
Fc_list = Fc_list*ScaleFactor;

%% Plotting
PPnames = {'F_y','theta_p','theta_pc','combined'};
PPcols = [1, 4, 5; 1, 6, 7; 1, 2, 3; 1, 8, 9];
for i = 1:4
    x1 = 100*[0, dy_list(PPcols(i,1)), dc_list(PPcols(i,1)), du_list(PPcols(i,1))];
    x1p5 = 100*[0, dy_list(PPcols(i,2)), dc_list(PPcols(i,2)), du_list(PPcols(i,2))];
    x2 = 100*[0, dy_list(PPcols(i,3)), dc_list(PPcols(i,3)), du_list(PPcols(i,3))];
    y1 = [0, Fy_list(PPcols(i,1)), Fc_list(PPcols(i,1)), 0];
    y1p5 = [0, Fy_list(PPcols(i,2)), Fc_list(PPcols(i,2)), 0];
    y2 = [0, Fy_list(PPcols(i,3)), Fc_list(PPcols(i,3)), 0];
    
    % plotting
    gcf = figure('Color',[1 1 1]);
    set(gcf,'units','inches','position',[1 1 1.5 1.5],'PaperUnits', 'Inches');
    gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
        'FontSize',10,...
        'FontName','Arial',...
        'Linewidth', 1,...
        'TickLength', [0.02 0.035],...
        'XLim',[0,10],'YLim',[0, 0.2]);%,...
    % 'XTick',[x_lim(1)+0.5,mean(x_lim),x_lim(2)-0.5]);%,...
    % 'YTick',[y_lim(1)+0.5,mean(y_lim),y_lim(end)-0.5]);
    set(gca, 'units','inches','position',[0.45 0.45 0.9 0.9])
    box(gca,'on');
    hold(gca,'all');
    plot(x1,y1,Lstyle,'Color',colors(i,:),'LineWidth',0.5)
    plot(x1p5,y1p5,Lstyle,'Color',colors(i,:),'LineWidth',1)
    plot(x2,y2,Lstyle,'Color',colors(i,:),'LineWidth',1.5)
    
    xlabel('Drift [%]','FontSize',10,'FontName','Arial')
    ylabel('Base Shear/Weight','FontSize',10,'FontName','Arial')
    print('-depsc', [pwd '/Figs/eps/backbone','_',struc,...
        '_',PPnames{i},'.eps']);
    if SAVE_FIG==true
        savefig([pwd '/Figs/eps/backbone','_',struc,...
            '_',PPnames{i},'.fig']);
    end
end
