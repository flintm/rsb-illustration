function M3_Plot_Hyperpolygon(gs,fX,Rzz,is,js,u_stars,alphas,Beta,name,ALPHA,STAR,TAN,SYS)
% gs = limit state functions
% fX = real space distributions
% Rzz = Nataf correlation coefficients
% is = variable X numbers (3 total)
% js = limit state function numbers (3 total)
% name is the alternative design configuration that generates this U-space
% and design point
%% plotting controls
c1 = [166,206,227]/255; % l blue
c2 = [251,154,153]/255; % pink
c3 = [51,160,44]/255; % green
c4 = [31,120,180]/255; % blue
c5 = [178,223,138]/255; % l green

% if in global mode, different colors
cz1 = [27,158,119]/255; % teal
cz2 = [217,95,2]/255; % orange
cz3 = [117,112,179]/255; % purple
%231,41,138
%102,166,30

cf = [0.2 0.2 0.20];
c = [c1; c2; c3; c4; c5];
cz = [cz1; cz2; cz3];
% #a6cee3 #1f78b4 #b2df8a #33a02c
alphF_f = 0.07;
alphF_h = 0.5;
alphE_h = 0.5;

view_ang = [-29.3800087169416 44.2945102064871];%[-41.94148245461929,26.502003001876194]; % [-1 -1 -1]; % 
u_plot_lim = 5;

pdf_vals = [1, 3, 5];
%% variables for transforming
u_lim = 15;
fz = makedist('Normal',0,1);
L = chol(Rzz)';
Linv = inv(L);
n = size(Rzz,1);

% for limit state functions
u = ones(n,1)*(-u_lim:0.5:u_lim);
[U1, U2, U3] = meshgrid(u(is(1),:),u(is(2),:),u(is(3),:));
xp = zeros(size(u));
zp = zeros(size(u));
up = zeros(size(u));
for i=1:n
   xp_min = icdf(fX{i},0.00001);
   xp_max = icdf(fX{i},1-0.00001);
   xp(i,:) = linspace(xp_min,xp_max,size(u,2));
end

for ii=1:size(u,2)
   z(:,ii) = L*u(:,ii); 
    for i=1:n
        xx(i,ii) = icdf(fX{i}, cdf(fz,z(i,ii)));
        zp(i,ii) = icdf(fz, cdf(fX{i},xp(i,ii)));
    end
    up(:,ii) = Linv*zp(:,ii);
end

[Z1, Z2, Z3] = meshgrid(z(is(1),:),z(is(2),:),z(is(3),:));
disp('Range of Z variables:')
disp([min(z(is(1),:)),max(z(is(1),:))])
disp([min(z(is(2),:)),max(z(is(2),:))])
disp([min(z(is(3),:)),max(z(is(3),:))])
[X1, X2, X3] = meshgrid(xx(is(1),:),xx(is(2),:),xx(is(3),:));
[XP1, XP2, XP3] = meshgrid(xp(is(1),:),xp(is(2),:),xp(is(3),:));
[UP1, UP2, UP3] = meshgrid(up(is(1),:),up(is(2),:),up(is(3),:));

% for contours
[ZZ1, ZZ2, ZZ3]  = meshgrid(-(u_plot_lim+1):0.1:(u_plot_lim+1),-(u_plot_lim+1):0.1:(u_plot_lim+1),-(u_plot_lim+1):0.1:(u_plot_lim+1));
[UU1, UU2, UU3] = meshgrid(-(u_plot_lim+1):0.1:(u_plot_lim+1),-(u_plot_lim+1):0.1:(u_plot_lim+1),-(u_plot_lim+1):0.1:(u_plot_lim+1));
fU = normpdf(UU1).*normpdf(UU2).*normpdf(UU3);
fZ = mvnpdf([ZZ1(:),ZZ2(:),ZZ3(:)], zeros(1,3), Rzz(is,is));
fZ = reshape(fZ,size(ZZ1));
%% make Z-space 3D contour figure
fCZ = figure('Color',[1,1,1],'Units','centimeters','Position',[1,1,15,15]);
axes('FontSize',12,'FontName','Arial','FontSmoothing','off',...
    'Xtick',[-u_lim:u_lim],'YTick',[-u_lim:u_lim],'Ztick',[-u_lim:u_lim],... 
    'TickDir','both',...%'XTickLabel','','YTickLabel','','ZTickLabel','',...
    'Xcolor','black','Ycolor','black','Zcolor','black',...
    'Linewidth',1); % 'Color','none',... 

hold all
hAxis = gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
daspect([1,1,1])
view(view_ang); 


for i=1:size(pdf_vals,2)
    isx = sqrt(pdf_vals(i)^2/3);
    isovalues(i) = mvnpdf(isx*ones(3,1));
    surf1 = isosurface(ZZ1,ZZ2,ZZ3,fZ,isovalues(i));
    p1 = patch(surf1);
    isonormals(ZZ1,ZZ2,ZZ3,fZ,p1);
    set(p1,'FaceColor',cf,'EdgeColor','none','FaceAlpha',alphF_f,'Clipping','off'); % set the color, mesh and transparency level of the surface
end

xlabel(['Z_', num2str(is(1))]);
ylabel(['Z_', num2str(is(2))]);
zlabel(['Z_', num2str(is(3))]);
xlim([-u_plot_lim,u_plot_lim]);
ylim([-u_plot_lim,u_plot_lim]);
zlim([-u_plot_lim,u_plot_lim]);
%% make U-space 3D contour figure
fC = figure('Color',[1,1,1],'Units','centimeters','Position',[1,1,15,15]);
axes('FontSize',12,'FontName','Arial','FontSmoothing','off',...
    'Xtick',[-u_lim:u_lim],'YTick',[-u_lim:u_lim],'Ztick',[-u_lim:u_lim],... 
    'TickDir','both',...%'XTickLabel','','YTickLabel','','ZTickLabel','',...
    'Xcolor','black','Ycolor','black','Zcolor','black',...
    'Linewidth',1); % 'Color','none',... 

hold all
hAxis = gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
daspect([1,1,1])
view(view_ang); 


for i=1:size(pdf_vals,2)
    isx = sqrt(pdf_vals(i)^2/3);
    isovalues(i) = mvnpdf(isx*ones(3,1));
    surf1 = isosurface(UU1,UU2,UU3,fU,isovalues(i));
    p1 = patch(surf1);
    isonormals(UU1,UU2,UU3,fU,p1);
    set(p1,'FaceColor',cf,'EdgeColor','none','FaceAlpha',alphF_f,'Clipping','off'); % set the color, mesh and transparency level of the surface
end

xlabel(['U_', num2str(is(1))]);
ylabel(['U_', num2str(is(2))]);
zlabel(['U_', num2str(is(3))]);
xlim([-u_plot_lim,u_plot_lim]);
ylim([-u_plot_lim,u_plot_lim]);
zlim([-u_plot_lim,u_plot_lim]);
%% set up u_star and alpha to plot
if STAR
    fS = figure;
    hold all
    for j=1:max(size(js))
        u_star = u_stars(js(j),is);
        plot3(u_star(1),u_star(2),u_star(3),'o','LineWidth',2,'Color', c(j,:))
    end
    fSZ = figure;
    hold all
    for j=1:max(size(js))
        z_star = L(is,is)*u_stars(js(j),is)';
        plot3(z_star(1),z_star(2),z_star(3),'o','LineWidth',2,'Color', c(j,:))
    end
end

if ALPHA
    fA = figure;
    hold all
    for j=1:max(size(js))
        u_star = u_stars(js(j),is);
        alpha = alphas(js(j),is);
        a = [u_star; u_star + alpha];
        plot3(a(:,1),a(:,2),a(:,3),'-','LineWidth',2,'Color', c(j,:))
    end
    fAZ = figure;
    hold all
    for j=1:max(size(js))
        z_star = transpose(L(is,is)*u_stars(js(j),is)');
        alpha = transpose(L(is,is)*alphas(js(j),is)');
        a = [z_star; z_star + alpha];
        plot3(a(:,1),a(:,2),a(:,3),'-','LineWidth',2,'Color', c(j,:))
    end
end
%% Original limit state surfaces Z space
fNZ = copyfig(fCZ);
disp('Making original limit state surfaces in Z-space')
for j=1:max(size(js))
    disp(num2str(js(j)))
    mesh = zeros(n,size(X1(:),1));
    mesh(is(1),:) = X1(:)';
    mesh(is(2),:) = X2(:)';
    mesh(is(3),:) = X3(:)';
    h = gs{js(j)}(mesh);
    h(isinf(h)) = NaN;
    if(all(h(~isnan(h))>0))
        disp(['   Only positive evaluations of g for LS',num2str(js(j))])
    end
    if(all(h(~isnan(h))<0))
        disp(['   Only negative evaluations of g for LS',num2str(js(j))])
    end
    h = reshape(h,size(X1));
    LSS = isosurface(Z1,Z2,Z3,h,0);
    pg = patch(LSS,'FaceColor', c(j,:),'EdgeColor','none',...
        'EdgeAlpha',alphE_h ,'FaceAlpha',alphF_h);
    isonormals(Z1,Z2,Z3,h,pg);
end

if ALPHA
    L = findobj(fAZ,'type','line');
    copyobj(L,findobj(fNZ,'type','axes'));
end

if STAR
    L = findobj(fSZ,'type','line');
    copyobj(L,findobj(fNZ,'type','axes'));
end

savefig([pwd,'/Figs/fig/Z_space_nonlinear_',...
    name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.fig']);
print('-depsc', [pwd,'/Figs/eps/Z_space_nonlinear',...
     name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.eps']);
%% Original limit state surfaces U space
fN = copyfig(fC);
disp('Making original limit state surfaces in U-space')
for j=1:max(size(js))
    disp(num2str(js(j)))
    mesh = zeros(n,size(X1(:),1));
    mesh(is(1),:) = XP1(:)';
    mesh(is(2),:) = XP2(:)';
    mesh(is(3),:) = XP3(:)';
    h = gs{js(j)}(mesh);
    h(isinf(h)) = NaN;
    if(all(h(~isnan(h))>0))
        disp(['   Only positive evaluations of g for LS',num2str(js(j))])
    end
    if(all(h(~isnan(h))<0))
        disp(['   Only negative evaluations of g for LS',num2str(js(j))])
    end
    h = reshape(h,size(X1));
    LSS = isosurface(UP1,UP2,UP3,h,0);
    pg = patch(LSS,'FaceColor', c(j,:),'EdgeColor','none',...
        'EdgeAlpha',alphE_h ,'FaceAlpha',alphF_h);
    isonormals(UP1,UP2,UP3,h,pg);
end

if ALPHA
    L = findobj(fA,'type','line');
    copyobj(L,findobj(fN,'type','axes'));
end

if STAR
    L = findobj(fS,'type','line');
    copyobj(L,findobj(fN,'type','axes'));
end

savefig([pwd,'/Figs/fig/U_space_nonlinear_',...
    name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.fig']);
print('-depsc', [pwd,'/Figs/eps/U_nonlinear_space_',...
     name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.eps']);
%% tangent planes
if TAN
    fT = copyfig(fC);
    disp('Making tangent plane surfaces')
    for j = 1:max(size(js))
        disp(num2str(j))
        [uty, utz] = meshgrid(u(1,:),u(3,:));
        t = 1/alphas(js(j),is(1))*(-alphas(js(j),is(2))*uty-alphas(js(j),is(3))*utz+Beta(js(j)));
        surf(t,uty,utz,'FaceColor',c(j,:),'EdgeColor','none',...
            'FaceAlpha',alphF_h,'EdgeAlpha',alphE_h,'Clipping','on')
    end
    
      
    if ALPHA
        L = findobj(fA,'type','line');
        copyobj(L,findobj(fT,'type','axes'));
    end
    
    if STAR
        L = findobj(fS,'type','line');
        copyobj(L,findobj(fT,'type','axes'));
    end
end
%% Systems 
if SYS
    disp('Making correlated Z-space surfaces')
    R = eye(3);
    R(1,2) = alphas(js(1),is)*alphas(js(2),is)';
    R(1,3) = alphas(js(1),is)*alphas(js(3),is)';
    R(2,3) = alphas(js(2),is)*alphas(js(3),is)';
    R(2,1) = R(1,2);
    R(3,1) = R(1,3);
    R(3,2) = R(2,3);
    
    z1 = (-u_lim:0.1:u_lim)';
    z2 = z1;
    z3 = z1;
    
    [ZZ1, ZZ2, ZZ3] = meshgrid(z1,z2,z3);
    ZZ = zeros([size(ZZ1),3]);
    ZZ(:,:,:,1) = ZZ1;
    ZZ(:,:,:,2) = ZZ2;
    ZZ(:,:,:,3) = ZZ3;
    fz = zeros(size(ZZ1));
    for i=1:size(ZZ,1)
        for j=1:size(ZZ,2)
            for k=1:size(ZZ,3)
                fz(i,j,k) = mvnpdf(reshape(ZZ(i,j,k,:),3,1),zeros(3,1),R);
            end
        end
    end
    
    figure('Color',[1,1,1],'Units','centimeters','Position',[1,1,15,15])
    axes('FontSize',12,'FontName','Arial','FontSmoothing','off',...
        'Xtick',[-u_lim:u_lim],'YTick',[-u_lim:u_lim],'Ztick',[-u_lim:u_lim],...
        'TickDir','both',...%'XTickLabel','','YTickLabel','','ZTickLabel','',...
        'Xcolor','black','Ycolor','black','Zcolor','black',... %     'Color','none',...
        'Linewidth',1,'clipping','off')
    hold all
    hAxis = gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
    hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
    hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
    daspect([1,1,1])
    view(view_ang);
    %camlight; lighting gouraud
    
    pdf_vals = [1, 3, 5];
    for i=1:size(pdf_vals,2)
        isx = sqrt(pdf_vals(i)^2/3);
        isovalues(i) = mvnpdf(isx*ones(3,1));
        %isovalues(i) = mvnpdf(pdf_vals(i)*ones(3,1),zeros(3,1),R);
        surf1 = isosurface(ZZ1,ZZ2,ZZ3,fz,isovalues(i));
        p1 = patch(surf1);
        isonormals(ZZ1,ZZ2,ZZ3,fz,p1);
        set(p1,'FaceColor',cf,'EdgeColor','none','FaceAlpha',alphF_f); % set the color, mesh and transparency level of the surface
    end
    u(:,u(1,:)<-u_plot_lim) = [];
    u(:,u(1,:)>u_plot_lim) = [];
    [uty, utz] = meshgrid(u(2,:),u(3,:));
    % use -Beta because of parallel system of cut sets--all must be
    % violated to cause system failure
    surf(-Beta(js(1))*ones(size(uty)),uty,utz,'FaceColor',cz1,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',alphE_h)
    [utx, utz] = meshgrid(u(1,:),u(3,:));
    surf(utx,-Beta(js(2))*ones(size(utx)),utz,'FaceColor',cz2,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',alphE_h)
    [utx, uty] = meshgrid(u(1,:),u(2,:));
    surf(utx,uty,-Beta(js(3))*ones(size(utx)),'FaceColor',cz3,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',alphE_h)
    
    xlabel(['Z_', num2str(js(1))]);
    ylabel(['Z_', num2str(js(2))]);
    zlabel(['Z_', num2str(js(3))]);
    xlim([-u_plot_lim,u_plot_lim])
    ylim([-u_plot_lim,u_plot_lim])
    zlim([-u_plot_lim,u_plot_lim])
    % export
    savefig([pwd,'/Figs/fig/Sys_Z_space_',name ,...
         name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.fig']);
    print('-depsc', [pwd,'/Figs/eps/Sys_Z_space_',name ,...
         name,'_',num2str(js(1)),'_',num2str(js(2)),'_',num2str(js(2)),'.eps']);
end
end