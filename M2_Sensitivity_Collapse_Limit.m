% Sensitivity of EAL/EL to collapse drift limit for Flint et al. 2019
% Written by Madeleine Flint, 2019-09-18.

close all
clc

%% Run main analysis
Sas_CS         = zeros(numGM,numStruct);
EDPs_CS        = zeros(numGM,numCase,numStruct);
pd_Cs_CS       = cell(numCase,numStruct,numFac);
paramNCs_CS    = zeros(3,numCase,numStruct,numFac);
fragSs_CS      = zeros(numSa,numDS,numCase,numStruct,numFac);
DamageProbs_CS = zeros(numSa,numDS+1,numCase,numStruct,numFac);
EL_PPs_CS      = zeros(numPPval,numPP,numStruct,numFac);
EL_PPs_cases_CS= zeros(numCase,numStruct,numFac);
disp('Performing PBEE assessment for all configurations of structure conditional on collapse drift limit:');
for l = 1:numStruct
    struc = strucs{l};
    disp(struc);
    for m = 1:numFac
        disp(['     ',num2str(m)]);
        RD_c  = RD_c_fac(m)*mu_RD_c(l);
        %importing Opensees results saved as matlab variables
        Sa = importdata(['Data/' IMname struc '.mat']);
        EDP = importdata(['Data/' EDPname struc '.mat']);
        
        [Sas_CS(:,l), EDPs_CS(:,:,l), pd, paramNCs_CS(:,:,l,m), ...
            fragSs_CS(:,:,:,l,m), DamageProbs_CS(:,:,:,l,m), ...
            EL_PPs_CS(:,:,l,m), EL_PPs_cases_CS(:,l,m)] = M2_PBEE_Simple(numCase, RD_c, ...
            hazardInfo, SaCalc, EDP, Sa, DSS(l,:),...
            DSNS(l,:), Repair_cost_S, Repair_cost_NS,...
            Repair_cost_C, T1(l), Years);
%         [Sas_CS(:,l), EDPs_CS(:,:,l),pd, ...
%             paramNCs_CS(:,:,l,m), fragSs_CS(:,:,:,l,m),DamageProbs_CS(:,:,:,l,m), ...
%             SaCalc, EL_PPs_CS(:,:,l,m), EL_PPs_cases_CS(:,l,m)] = ...
%                 M2_PBEE_Simple(struc, RD_c, ... 
%                 hazname, EDPname, IMname, Years);
         pd_Cs_CS(:,l,m) = pd;
    end
end

%% Caculate distribution parameters for X3, P(Collapse|Sa_MCE)
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.MaxFunEvals = 10000;
opts.MaxIter = 8000;
opts.Robust = 'LAR';
P_c_RD_c = zeros(numCase,numStruct,numFac);
P_c_fit  = cell(numCase,numStruct);
gof      = cell(numCase,numStruct);
P_c      = zeros(numCase,numStruct);
ELLR_c   = zeros(numCase,numStruct);
EL_PPs_c = zeros(numPPval,numPP,numStruct);

disp('Fitting conditional collapse distributions');
for l = 1:numStruct
    struc = strucs{l};
    disp(struc);
    pd_RD_c  = makedist('lognormal','mu',log(mu_RD_c(l)),...
        'sigma',sigma_RD_c(l));
    rd = (mu_RD_c(l)*RD_c_fac)';
    RD = [0.001; rd; 2*rd(numFac)];
    RD_mid = 0.5*RD(1:(numFac+1))+0.5*RD(2:(numFac+2));
    pr_RD_c = cdf(pd_RD_c, RD_mid(2:(numFac+1)))-...
        cdf(pd_RD_c, RD_mid(1:numFac));
    for k = 1:numCase
        disp(['   case: ' num2str(k)]);
        for m= 1:numFac
                P_c_RD_c(k,l,m)  = cdf(pd_Cs_CS{k,l,m},log(Sa_MCE(l)));
        end
        [P_c_fit{k,l}, gof{k,l}] = fit(rd, reshape(P_c_RD_c(k,l,:),numFac,1), ft, opts );
        P_c(k,l) = pr_RD_c'*reshape(P_c_RD_c(k,l,:),numFac,1);
        ELLR_c(k,l) = pr_RD_c'*reshape(EL_PPs_cases_CS(k,l,:),numFac,1);
    end
    EL_PPs_c(:,:,l) = [ELLR_c(1,l)*ones(1,numPP);
                       ELLR_c(2:2:numCase,l)';
                        ELLR_c(3:2:numCase,l)'];
end

g3_formula = subs(str2sym([num2str(a(3,3)),'*(',...
    formula(P_c_fit{1,1}),')+',num2str(a(3,8))]),'x','x3');
grad_g3_formula = diff(g3_formula,1);
%% Evaluate correlation of collapse drift limit and drift at DBE, MCE
% This is very computationally expensive and not advised
if COMPUTE_CORR_MC
    close all
    rho_x1_x3 = zeros(numCase,numStruct);
    rho_x2_x3 = zeros(numCase,numStruct);
    disp('computing correlation between MCE drift and collapse limit');
    n_rnd = 1000000;
    rand_u1 = rand([n_rnd,1]);
    disp('limit side');
    for l = 1:numStruct
        disp(strucs{l})
        pd_RD_c  = makedist('lognormal','mu',log(mu_RD_c(l)),...
            'sigma',sigma_RD_c(l));
        d_lim(:,l) = icdf(pd_RD_c,rand_u1);
        d_lim_sort(:,l) = sort(d_lim(:,l));
    end
    disp('demand side; cases');
    rand_u2 = rand([n_rnd,1]);
    for k=5:numCase
        disp(k);
        figure
        hold on
        for l = 1:numStruct
            disp(['     ' strucs{l}]);
            C = NaN*ones(n_rnd,1);
            %d_DBE = icdf(f_RD_Sa_DBE{k,l},rand_u3);%  no collapses ever found
            d_MCE = icdf(f_RD_Sa_MCE{k,l},rand_u2);
            d_MCE_sort = sort(d_MCE);
            done = false;
            i1 = 1;
            in = n_rnd;
            ip5 = floor((i1+in)/2);
            is  = [i1; ip5; in]
            it = 0;
            while ~done
                temp = transpose(d_MCE_sort)>=d_lim_sort(is,l);
                C(is)  = sum(temp,2);
                [row1, col1] = find(temp,1);
                [row2, col2] = find(temp,1,'last');
                if (any(size(row1)==0)) % all smaller
                    C(i1:in) = 0;
                else
                    if row2 == 3 % need to fill in all values up to in
                        C(i1) = n_rnd - col1 + 1;
                        C(ip5) = sum(temp(2,:));
                        C(in) = sum(temp(3,:));
                        for i = (i1+1):(in-1)
                            C(i) = sum(transpose(d_MCE_sort(col1:col2))>=d_lim_sort(i,l));
                        end
                    else % know that row2==1 or row2==2
                        [mx, wh] = max(C(is)==0); 
                        C(is(max(1,wh-1)):in) = 0;
                    end
                end
                i1 = max(1,find(isnan(C),1));
                if any(size(i1)==0)
                    done=true;
                else
                    in = find(isnan(C),1,'last');
                    ip5 = floor((i1+in)/2);
                    is = [i1; ip5; in]
                end
                it = it + 1;
                if it>100
                    error('Finished by iteration limit');
                end
            end
            
            % compute correlation / fit model
            P_c_k = C/n_rnd;
            tb = table(d_lim_sort(:,l),P_c_k,'VariableNames',{'DriftLim','ProbColl'});
            c_fit = fitlm(tb,'ProbColl~DriftLim');
            rhos = corr([d_lim_sort(:,l),P_c_k]); %C_d,
            rho_x2_x3(k,l) = rhos(1,2);
            rho_x1_x3(k,l) = 0;
            
            % make plots
            subplot(2,2,2*(l-1)+1)
            plot(d_MCE, d_lim(:,l),'.k');
            xlabel('Drift at MCE');
            ylabel('Collapse drift limit');
            hold on;
            plot([0;1],[0;1],'-');
            xlim([0,.08]);
            ylim([0,.08])
            axis square;
            subplot(2,2,2*l)
            plot(c_fit);
            hold on;
            box on;
            if max(P_c_k)>0
                ylim2 = 10^(ceil(log10(max(P_c_k))));
            else
                ylim2 = 1e-05;
            end
            xlim([0,.08]);
            ylim([0,ylim2]);
            text(0.04, ylim2/4, ['Rho_x2_x3 = ' sprintf('%0.2f',rho_x2_x3(k,l))]);
        end
        sgtitle(['Case ' num2str(k)]);
        if SAVE_FIG==true
            savefig([pwd,'/Figs/fig/Drift_Collapse_corr_',num2str(k),'v2.fig']);
        end
        print('-depsc', [pwd '/Figs/eps/Drift_Collapse_corr_',num2str(k),'v2.eps']);
    end
    for l=1:numStruct
        r_x2_x3(:,l) = round(rho_x2_x3(:,l),1);
        r_x2_x3(isnan(r_x2_x3(:,l)),l) = 0;
        r_x1_x3(1:numCase,l) = 0;
    end
    save('Data/drift_collapse_corr_v2.mat','rho_x2_x3', 'r_x2_x3', 'r_x1_x3') 
else
    for l=1:numStruct
        disp(strucs{l})
        for k=1:numCase
            disp(num2str(k))
            rd = (mu_RD_c(l)*RD_c_fac)';
            dat = reshape(P_c_RD_c(k,l,:),numFac,1);
            rho_x2_x3(k,l) = corr(rd,dat);
        end
        r_x2_x3(:,l) = round(rho_x2_x3(:,l),1);
        r_x1_x3(:,l) = 0;
    end
save('Data/drift_collapse_corr_v3.mat','rho_x2_x3', 'r_x2_x3', 'r_x1_x3') 
end
%% Smoothed collapse conditional distributions
l = 1; 
k = 1; %base
struc = strucs{l};
rd = (mu_RD_c(l)*RD_c_fac)';
dat = reshape(P_c_RD_c(k,l,:),numFac,1);
rd_smth = 0.03:0.001:0.1;
smth = P_c_fit{k,l}(rd_smth);
gcf = figure('Color',[1 1 1]);
set(gcf, 'units','inches','position',[1 1 3 2],'PaperUnits', 'Inches');
gca1 = axes('Parent',gcf,'YGrid','off','XGrid','off',...
    'FontSize',10,...
    'FontName','Arial',...
    'Linewidth', 1,...
    'TickLength', [0.02 0.035],...
    'YLim',[0, 0.2],'XLim',[0.03,0.1],'Xcolor',[0,0,0],'Ycolor',[0,0,0],...
    'XTick',[0.04,0.06,0.08,0.1]);%,...
    %'YTick',[0.1, 0.3, 0.5, 0.7],'YTickLabel',{'0.1%','0.3%','0.5%','0.7%'});
box(gca1,'on');
hold(gca1,'all');
hold all
Lstyle = Lstyles{l};
Mstyle = 'o';

plot(rd,dat,Mstyle,'Color',[0,0,0])
plot(rd_smth,smth,Lstyle,'Linewidth',1.5,'Color',[0,0,0]);%,'MarkerEdgeColor',...

xlabel('Collapse Drift Limit, RD_C');
ylabel('P(Collapse|Sa_M_C_E,RD_C)');
if SAVE_FIG
    savefig([pwd,'/Figs/fig/Collapse_fit_',struc,'.fig']);
end
print('-depsc', [pwd '/Figs/eps/Collapse_fit_',struc,num2str(k),'.eps']);
%% conditional EL values; concrete and steel separate
for l = 1:numStruct
    gcf(l) = figure('Color',[1 1 1]);
    set(gcf(l), 'units','inches','position',[1 1 3 2],'PaperUnits', 'Inches');
    struc = strucs{l};
    % plot EL as a function of drift limit for collapse
    gca1 = axes('Parent',gcf(l),'YGrid','off','XGrid','off',...
        'FontSize',10,...
        'FontName','Arial',...
        'Linewidth', 1,...
        'TickLength', [0.02 0.035],...
        'YLim',[0, 0.9],'XLim',[0.03,0.1],...
        'XTick',[0.04,0.06,0.08,0.1],...
         'YTick',[0.1, 0.3, 0.5, 0.7, 0.9],'YTickLabel',{'0.1%','0.3%','0.5%','0.7%','0.9%'});
    box(gca1,'on');
    hold(gca1,'all');
    hold all
    Lstyle = Lstyles{l};
    for kk = 1:numPP 
        Mstyle = MstylesComb{kk};
        for k = 1:numPPval 
            plot((mu_RD_c(l)*RD_c_fac)',reshape(EL_PPs_CS(k,kk,l,:),numFac,1),...
                Lstyle,'Linewidth',LwidthsL(k),'Color',colors(kk,:));%,'MarkerEdgeColor',...
               % colors(j,:), ,'MarkerSize',2,...
                %'Marker',Mstyle)
        end
    end
    xlabel('Collapse Drift Limit');
    ylabel('E[Lifetime Loss Ratio]');
    if SAVE_FIG
        savefig([pwd,'/Figs/fig/EL_Collapse_Limit_',struc,'.fig']);
    end
    print('-depsc', [pwd '/Figs/eps/EL_Collapse_Limit_',struc,'.eps']);
    
    % plot assumed drift limit distribution
    gcf2(l) = figure('Color',[1 1 1]);
    set(gcf2(l), 'units','inches','position',[1 1 3 1],'PaperUnits', 'Inches');
    gca2 = axes('Parent',gcf2(l),'YGrid','off','XGrid','off',...
        'Position',[0.197129637058036,0.20375,0.707870362941964,0.72125],...
        'FontSize',10,...
        'FontName','Arial',...
        'Linewidth', 1,...
        'TickLength', [0.02 0.035],...
        'XLim',[0.03,0.1],'YLim',[0, 80],...
        'XTickLabel','','XTick',[0.04,0.06,0.08,0.1],...
        'YTick',[0,20,40,60, 80],...
        'YDir','reverse');%,... %
         
    box(gca2,'on');
    hold(gca2,'all');
    RD = 0.03:0.001:0.1;
    f_RD = lognpdf(RD, log(mu_RD_c(l)), sigma_RD_c);
    plot(RD,f_RD,...
        Lstyle,'Linewidth',L_widths(k),'Color','black')
    
    ylabel('f_C')
    if SAVE_FIG
        savefig([pwd,'/Figs/fig/f_Collapse_Limit_',struc,'.fig']);
    end
    print('-depsc', [pwd '/Figs/eps/f_Collapse_Limit_',struc,'.eps']);
end
%% EL values with collapse uncertainty integrated through, conc+steel
gcf = figure('Color',[1 1 1]);
set(gcf, 'units','inches','position',[1 1 3 3],'PaperUnits', 'Inches');
gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
            'FontSize',10,...
    'FontName','Arial',...
    'Linewidth', 1,...
    'TickLength', [0.02 0.035],... 
    'YLim',[0, 0.9],...
    'XTick',[0,0.5,1],'Xcolor',[0,0,0],...
    'YTick',[0.1, 0.3, 0.5, 0.7,0.9],'Ycolor',[0,0,0],...
    'YTickLabel',{'0.1%','0.3%','0.5%','0.7%','0.9%'});
box(gca,'on');
hold(gca,'all');
MstylesComb = {'o';'s';'^';'V'};
for l = 1:numStruct
    Lstyle = Lstyles{l};
    for k = 1:numPP
        Mstyle = MstylesComb{k};
        plot(theta,EL_PPs_c(:,k,l,1),Lstyle,'MarkerEdgeColor',colors(k,:),...
            'Marker',Mstyle,...
             'Color',colors(k,:))
    end
end
xlabel('Performance Parameter');
ylabel('Expected Lifetime Loss Ratio (all RD_c)');
if SAVE_FIG
    savefig([pwd,'/Figs/fig/EL_conc_steel_RD_c_all.fig']);
end
print('-depsc', [pwd,'/Figs/eps/EL_conc_steel_RD_c_all.eps']);

