% Transforms a pushover curve to a trilinear SDOF backbone
% written by Mohsen Zaker-Esteghamati
clear
clc
pushdata=csvread('4story_concrete.csv');
W=(max(pushdata(:,2))/1.6)/0.092;
pushdata(:,2)=pushdata(:,2)./W;

%% 
xArea=min(pushdata(:,1)):0.001:max(pushdata(:,1));
yArea=interp1(pushdata(:,1),pushdata(:,2),xArea);
PushArea=trapz(xArea,yArea);
%% 
%to get the linear part slope
m=1;
k=(2*pi/1.16)^2*m;
Fy=max(pushdata(:,2));
xy=(Fy/k);
%defining the difference between curves 
fun=@(x)abs((0.5*Fy*xy+0.5*(1+x(1))*Fy*(x(2)-(xy))+0.5*(pushdata(end,1)-x(2))*Fy*x(1))-PushArea);
%initial guess
X0=[0.9,0.01];
%setting upperbond and lower bond for the solution
[param,fval,~,~,~,~,~]=fmincon(fun,X0);
% Plotting the results
%% 

Xp=[0, xy, param(2), 0.1];
Yp=[0,Fy, Fy*param(1), 0];
xscaled_y_ind=find(pushdata(:,2)<0.465*Fy&pushdata(:,2)>0.45*Fy);
xscaled=Xp.*[0.46,0.46,0.46,0.46];
yscaled=Yp.*0.46;
% 
% xscaled=[0 pushdata(xscaled_y_ind,1) param(2) 0.1];
% yscaled=[1 0.46*Fy 0.46*Fy*param(1) 0];
%% 
% Xp=[0	0.004501518	0.044973019	0.089347028]
% Yp=[0	0.277005556	0.308884318	0]
% % pushdata(:,1)=pushdata(:,1)./(15+3*13);
% % pushdata(:,2)=pushdata(:,2)./(193/0.092);
box on
set(gcf,'units','inches','position',[1 1 2 2],'PaperUnits', 'Inches', 'PaperSize', [2, 2])
plot(pushdata(:,1),pushdata(:,2),'k','lineWidth',1.2)
hold on
plot(Xp,Yp,'--r','lineWidth',1.5)
plot(xscaled,yscaled,'lineWidth',1.5,'color',[0.6 0.6 0.6])

set(gca,'FontSize',10)
set(gca,'FontName','Times')
xlabel('Roof drift')
ylabel('V/W')
ylim([0 0.24])
yticks([0:0.12:24])
xlim([0 0.1])
% lg=legend('Pushover curve','SDOF Backbone fit')
% currentLegendPosition = lg.Position;
% newLegendPosition = [0.5 0.7 currentLegendPosition([3 4])];
% % Set new position
% lg.Position = newLegendPosition;
% legend boxoff  
%% scaling
Xp=[0, param(1,1), param(2,1), param(3,1)]./(15+3*13);
% Yp=[0,param(1,2), param(2,2), param(3,2)]./(193/0.092);

        
