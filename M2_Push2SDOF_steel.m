%this script trasnform a pushover curve to a trilinear SDOF backbone
clear
clc
pushdata=csvread('4story_steel.csv');

%% 
i=1;
while i<length(pushdata(:,1))
    if pushdata(i,1)==pushdata(i+1,1)
        pushdata(i,:)=[];

    end
    i=i+1;
end

xArea=min(pushdata(:,1)):0.001:max(pushdata(:,1));
yArea=interp1(pushdata(:,1),pushdata(:,2),xArea);
PushArea=trapz(xArea,yArea);
%% 

%to get the linear part slop
k=(2*pi/1.32)^2;
Fy=max(pushdata(:,2));

%defining the difference between curves 
fun=@(x)abs((0.5*Fy*(Fy/k)+0.5*(1+x(1))*Fy*(x(2)-(Fy/k))+0.5*(pushdata(end,1)-x(2))*Fy*x(1))-PushArea);
% fun=@(X)abs(PushArea-((X(1,1)*a*X(1,1))*0.5+(a*X(1,1)+1.05*X(1,2))*0.5*(abs(X(2,1)-X(1,1)))+0.5*(abs(ultXpush-X(2,1)))*1.05*X(1,2)));
%initial guess
X0=[0.95,0.018];
%setting upperbond and lower bond for the solution


[param,fval]=fmincon(fun,X0);
% Plotting the results
%% 
x_y=(Fy/k);
Xp=[0, x_y, param(2), pushdata(end,1)];
Yp=[0,Fy, Fy*param(1), 0];
xscaled=Xp.*[0.46,0.46,0.46,0.46];
yscaled=Yp.*0.46;
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

        
