% Main script to run analyses and create figures in support of Flint et al. 
% 2019, "A PERFORMANCE-BASED DECISION SUPPORT SYSTEM FOR EARLY DESIGN OF 
% MULTI-HAZARD RESILIENT, SUSTAINABLE BUILDINGS"
% Written by Madeleine Flint, 2019-09-18 to 2019-11-13
% Requires that Data and Fig subfolders are available
% Suggest using Matlab "Run and Advance" option chunk-by-chunk
clear
close all
clc
SAVE_FIG = false; % will save .fig in addition to .eps

%% 1.  Setup: define structures, performance parameters, damage states
% Two structure type options
strucs = {'conc', 'steel'};
numStruct = size(strucs,2);

% structural properties:
% Ratio of S1 Charleston vs. design location
% concrete=0.344/0.742=0.4636; for steel=0.304/0.656=0.4634
T1 = [1.16; % period of concrete structure
      1.32]; % period of steel structure
backbone_SF = [0.46; 0.46]; 
mu_RD_c     = 1.5*[0.0370; 0.0471]; % multiple of Hazus high-code 'complete'
sigma_RD_c  = [0.1; 0.1]; % dispersion
Sa_DBE      = backbone_SF.*[0.742; 0.656]; % Sa_DBE original from USGS
Sa_MCE      = backbone_SF.*[0.984; 1.113];  % Sa_MCE original from USGS
numFac      = 50;
RD_c_fac    = 0.5:(1/(numFac-1)):1.5; % for collapse uncertainty
numFac      = size(RD_c_fac,2);

% impacts 
c_init    = [8.64, 8.91];    % Cost, Millions USD
e_init    = [28.4, 18.2];    % Embodied energy, TJ
e_op(1,:) = 606*[1, 1.02, 0.96]; % Operational energy, TJ, theta 0.5, 0, 1
e_op(2,:) = 597*[1, 1.015, 0.96];  
Years     = 50;           % Building lifetime
c_energy  = 0.0028;    % Energy cost, M$/TJ
ELMR      = 0.56; % Expected lifetime maintenance cost ratio for 50 years
CV        = [NaN; NaN; NaN; 0.05; NaN; 0.05; 0.2]; % assumed coeff. of var.
n_x       = size(CV,1);

% Embodied energy, concrete TJ
%  foundation 4.98
%  structure  7.81
%  envelope   3.87
%  floors     6.65
%  roofs      2.6
%  partions   2.50
% 28.41 TJ total

% Embodied energy, steel, TJ
%  foundation 4.98
%  structure  4.31
%  envelope   3.87
%  floors     2.53
%  roof       0.04
%  partions   2.50
%  18.23 TJ total
e_sfleh = [0 0; 
			4.98 4.98; 
			7.81 4.31;
		    3.87 3.87; 
			6.65+2.6, 2.53+.04]./e_init; % % of total
e_ns = (3.87 + 2.5)./e_init;

% Lateral system performance parameters and values & cost escalation
% increased cost/embodied energy also provided for envelope PPs related
% to window-wall ratio, although those cases are considered only at the end
% of the script
PP = transpose([1, 1.5, 2]);
theta = (PP-min(PP))/(max(PP)-min(PP));
numPPval = size(PP,1);
PPnames = {'F_y','theta_p','theta_pc','combined'};
numPP  = size(PPnames,2);
PPcols = [1, 4, 5; 1, 6, 7; 1, 2, 3; 1, 8, 9]; % for PBE_Simple
numCase = (numPPval-1)*numPP+1;

caseNames   = cell(numCase,1);
configNames = cell(numCase,numStruct);
for l=1:numStruct
    struc = strucs{l};
    configNames{1,l} = [struc,', Base'];
    caseNames{1} = '--';
   for kk=1:numPP
       for k=2:numPPval
           caseNames{2*kk+k-2} = [PPnames{kk},': ',...
                num2str(theta(k))];
            configNames{2*kk+k-2,l} = [struc,', ',PPnames{kk},': ',...
                num2str(theta(k))];
       end
   end
end

% cost and energy escalation to improve performance
% rows are PP, columns are PP val, 3rd dim is structure type
% Fy scales with performance increase on lateral system cost/energy
% theta_p scales at 50% on lateral
% theta_pc scales at 5% on gravity
c_PP_l(:,:,1) = [0.15 0.30; 0.08, 0.15; 0.04, 0.07]; % concrete, M$
c_PP_l(:,:,2) = [0.13 0.27; 0.07, 0.13; 0.03, 0.06]; % steel, M$
e_PP_l(:,:,1) = e_init(1)*[e_sfleh(3,1)*[0.5 1; 0.25 0.5];...
    e_sfleh(5,1)*[0.05 0.1]]; % concrete EE proportions, struc + floor TJ
e_PP_l(:,:,2) = e_init(2)*[e_sfleh(3,2)*[0.5 1; 0.25 0.5];...
   e_sfleh(5,2)*[0.05 0.1]]; % steel, struc + floor TJ

% combined for lateral
for l=1:numStruct
    c_case(:,l) = [0; ...
        reshape(transpose(c_PP_l(:,:,l)),(numPPval-1)*(numPP-1),1);...
        transpose(sum(c_PP_l(:,:,l),1))];
    e_case(:,l) = [0; ...
        reshape(transpose(e_PP_l(:,:,l)),(numPPval-1)*(numPP-1),1);...
        transpose(sum(e_PP_l(:,:,l),1))];
end

% Envelope window-wall ratio PP values & escalation
% base case corresponds to PP_e of 0.5, 31%, 0->40%, 1->22%
PP_e = [0.5 0 1]; 
PP_e_names{1} = ', PP_e: 0.5';
PP_e_names{2} = ', PP_e: 0';
PP_e_names{3} = ', PP_e: 1';
numPPeVal = size(PP_e,2);
c_PP_e = [0, 0; 0.048, 0.048; -0.048, -0.048]; % M$
%O/W Ratio	22%	27%	31%	35%	40% (PP 1 to 0, but using only 22, 31, 40)
%Concrete	23.514174	23.435817	23.373125	23.310479	23.232053
%steel	17.957833	17.924610	17.897991	17.871372	17.838099
e_PP_e(:,:) = [0 0;  -0.14 -0.07; 0.14 0.06]; %% TJ

% putting together full configuration names
names = cell(numCase,numStruct,numPPeVal);
namesNS = cell(numCase,numPPeVal);
for l = 1:numStruct
    for p = 1:numPPeVal
        for k = 1:numCase
            names{k,l,p} = [configNames{k,l},PP_e_names{p}];
            namesNS{k,p} = [caseNames{k},PP_e_names{p}];
        end
    end
end
names = names(:);
namesNS = namesNS(:);

% Information relevant to PBEE analysis / loading results
% Cloud Analysis results
numGM = 80;
SaCalc = (0:0.001:2)';
numSa = size(SaCalc,1); % hard-coded to match M2_PBEE_Simple
EDPname = 'EDP_';
IMname = 'Sa_';
hazname = 'hazardInfo';
% Damage states: structural and drift-sensitive non-structural from Hazus
% based on HAZUS high-code, using roof drifts
DSnames = {'slight','moderate','extensive','complete','collapse'};
numDS   = size(DSnames,2) - 1; % Hazus-MH of interest
DSS(1,:)  = [1.5,3,9,24]/648; % concrete, 648 = roof height, inches
DSS(2,:)  = [2.16,4.32,10.80,28.80]/612;% steel, 612 = height, inches
DSNS(1,:) = [1.80,3.60,11.25,22.50]/648; % concrete
DSNS(2,:) = [2.16, 4.32, 13.50, 27.0]/612; %steel
% as % of total replacement cost
Repair_cost_S  = [0.4,1.9,9.6,19.2]; % 4 structural
Repair_cost_NS =[0.7 ,3.3, 16.4, 32.9]; % 4 non-structural drift-sensitive
Repair_cost_C = 110;
r_x1_x2 = [0.8, 0.5; % correlation of DBE & MCE drift for cases, struct
           0.7, 0.7; % Fy
           0.8, 0.7;
           0.8, 0.8; % theta_p
           0.9, 0.8;
           0.9, 0.8; % theta_pc
           0.9, 0.8;
           0.9, 0.7; % comb
           0.9, 0.7];

% placeholder distributions
pds = cell(n_x,numStruct);
pds(:,1) = ...
      {'lognormal'; % RD at Sa_DBE
       'lognormal'; % RD at Sa_MCE
       'lognormal'; % Collapse drift limit
       'normal';    % Initial cost
       'normal';    % Life-cycle cost
       'normal';    % Embodied energy
       'normal'};   % Operational Energy
for i=2:numStruct
    pds(:,i) = pds(:,1); % use same distribution type for all structures
end
paramNames = cell(n_x,2,numStruct);
paramNames(:,1,:) = {'mu'};
paramNames(:,2,:) = {'sigma'};
paramVals = zeros(n_x,2,numCase,numStruct,numPPeVal);
%% 2.  Error encoding setup
R_adj = 0.001;  % nudge to correlation values if rank is deficient
R_adj_max = 60; % don't allow adjustments larger than R_adj*R_Adj*max
pC_min = 0.01;  % skip calculating pf3 if likelihood of collapse from 
                % convolution is below this value

record_Rzz_err = zeros(0,3);
record_g3_err = zeros(0,3);
record_R_UU_3_err = zeros(0,3);
record_R_UU_4_err = zeros(0,3);
%% 3.  Plotting setup
%defining color for plotting
color1 = '#42d483'; % theta_pc, green
color2 = '#f5c842'; % mustard, F_y
color3 = '#12a2db'; % blue, theta_p
color4 = '#25a88e'; % teal, combined
colors = [color2; color3; color1; color4];

% Linestyles and markers
Lstyles = {'-','--'}; % concrete; steel
Lwidths = [1, 0.5, 1.5]; % envelope PP value
Mstyles = {'o','s'}; % concrete; steel
MstylesComb = {'o';'s';'^';'V'}; % Fy, theta_p, theta_pc, combined
LwidthsL = [0.5, 1, 1.5]; % PP of 1, 1.5, 2

% x labels
xlabs = [true, true];
%% 4.  Module 0: Decision preferences & indifference curves
% Limit state functions (indifference curves)
% m x (n+1) vector of linear LSF coefficients
a =[-1,0,0,0,0,0,0,0.02;... % 1 DBE Drift
    0,-1,0,0,0,0,0,0.04;... % 2 MCE Drift
    0,0,-1,0,0,0,0,0.1;...  % 3 MCE Collapse probability
    0,0,0,-1,0,0,0,9;...  % 4 Initial Cost 
    0,0,0,0,-0.33,-0.13,0,8;... % 5 LCC versus embodied energy
    0,0,0,-1,0,0,-0.01,15.1;...  % 6 IC vs operation energy
    0,0,-5.5,-1,0,0,0,8.55;...   % 7 MCE collapse versus IC
    0,0,0,-1,-0.5,0,0,14.75;    % 8 IC versus LCC
    -1700,0,0,0,0,0,-0.01,20;...% 9 DBE Drift vs. Operational Energy
    0,-6000,0,0,0,0,-0.38,25];% 10 MCE Drift vs. Operational Energy
m_g  = size(a,1);      % number of limit state functions
g  = sym('g',[m_g,1]);
%grad_g = sym('grad_g',[m,1]);
gs = cell(m_g,1);
grad_gs = cell(m_g,1);

x   = sym('x',[n_x,1]); % Random variables

disp('Creating MATLAB functions for limit states:');
for j = [1:2,4:m_g] % g3 is non-linear
   disp(['Limit State ', num2str(j)]);
   g(j)  = a(j,1:n_x)*x + a(j,n_x+1);
   gs{j} = matlabFunction(g(j),'vars',{x});
   grad_g = sym(transpose(a(j,1:n_x)));
   grad_gs{j} = matlabFunction(grad_g,'vars',{x});
end
%% 5.  Module 2: Run main performance-based engineering analysis
Sas         = zeros(numGM,numStruct);
EDPs        = zeros(numGM,numCase,numStruct);
pd_Cs       = cell(numCase,numStruct);
paramNCs    = zeros(3,numCase,numStruct);
fragSs      = zeros(numSa,numDS,numCase,numStruct);
DamageProbs = zeros(numSa,numDS+1,numCase,numStruct);
EL_PPs      = zeros(numPPval,numPP,numStruct);
EL_PP_cases = zeros(numCase,numStruct);
EAL_NS      = zeros(numCase,numStruct);
EAL_S       = zeros(numCase,numStruct);
EAL_C       = zeros(numCase,numStruct);
f_RD_Sa_DBE = cell(numCase,numStruct);
f_RD_Sa_MCE = cell(numCase,numStruct);
% Main analysis to obtain expected loss
load(['Data/' hazname '.mat']); % variable must be hazardInfo
disp('Performing PBEE assessment for all configurations of structure:');
for l=1:numStruct
    struc = strucs{l};
    disp(struc);
    
    %importing Opensees results saved as matlab variables    
    Sa = importdata(['Data/' IMname struc '.mat']);
    EDP = importdata(['Data/' EDPname struc '.mat']);
    RD_c  = mu_RD_c(l);
    
    [Sas(:,l), EDPs(:,:,l), pd_Cs(:,l), paramNCs(:,:,l), ...
        fragSs(:,:,:,l), DamageProbs(:,:,:,l), ...
        EL_PPs(:,:,l), EL_PP_cases(:,l), EAL_S(:,l),...
        EAL_NS(:,l), EAL_C(:,l)] = M2_PBEE_Simple(numCase, RD_c, ...
                                  hazardInfo, SaCalc, EDP, Sa, DSS(l,:),...
                                  DSNS(l,:), Repair_cost_S, Repair_cost_NS,...
                                  Repair_cost_C, T1(l), Years);
    for k = 1:numCase
        f_RD_Sa_DBE{k,l} = makedist('lognormal',...
            'mu',paramNCs(1,k,l)+paramNCs(2,k,l)*log(Sa_DBE(l)),...
            'sigma',paramNCs(3,k,l));
        f_RD_Sa_MCE{k,l} = makedist('lognormal',...
            'mu',paramNCs(1,k,l)+paramNCs(2,k,l)*log(Sa_MCE(l)),...
            'sigma',paramNCs(3,k,l));
    end
end
%% 6.  Module 2: Collapse sensitivity and EDP/DS distributions.
% If the Monte Carlo analysis is to be re-run, the value of
% COMPUTE_CORR_MC at the top of M2_Sensitivity_Collapse_Limit.m
% must be changed to true.
COMPUTE_CORR = 'load'; % 'MC' 'load'
if strcmp(COMPUTE_CORR,'disc') || strcmp(COMPUTE_CORR,'MC')
    M2_Sensitivity_Collapse_Limit
    save('Data/drift_distributions.mat','f_RD_Sa_DBE','f_RD_Sa_MCE','P_c_fit',...
     'P_c','g3_formula','grad_g3_formula','ELLR_c','EL_PPs_c');
else
    if strcmp(COMPUTE_CORR,'load')
        load('Data/drift_distributions.mat')
        load('Data/drift_collapse_corr_disc.mat') % MC version also provided
    else
       error('Must compute correlations through Monte Carlo, discretized values, or load'); 
    end
end
%% 7.  Module 2: Calculate parameters for probabilistic distributions
% Random variables (decision attributes/metrics)
fX  = cell(n_x,numCase,numStruct,numPPeVal);
x0  = zeros(n_x,numCase,numStruct,numPPeVal);
Rzz = zeros(n_x,n_x,numCase,numStruct,numPPeVal);
R_err = 0;
R_err_max = 0;
for m = 1:numPPeVal
    for k=1:numCase
        for l = 1:numStruct
            % central values (means/medians)ç
            paramVals(:,1,k,l,m)  = ...
                [NaN;
                NaN;
                mu_RD_c(l);
                c_init(l) + c_case(k,l) + c_PP_e(m); % initial cost
                (c_init(l) + c_case(k,l) + c_PP_e(m))*...
                (1+ELMR+ELLR_c(k,l)/100)+c_energy*e_op(l,m);% LCC
                (e_init(l)+e_case(k,l) + e_PP_e(m))*...
                (1+ELLR_c(k,l)/100); % embodied E
                e_op(l,m)];          % operational energy
            
            % scale (standard deviations/dispersions)
            paramVals(:,2,k,l,m)  = [NaN;
                NaN;
                sigma_RD_c(l);
                CV(4)*paramVals(4,1,k,l,m); 
                NaN;
                CV(6)*paramVals(6,1,k,l,m); 
                CV(7)*paramVals(7,1,k,l,m)];
            paramVals(5,2,k,l,m) = (1+ELMR+EL_PP_cases(k,l)/100)*...
                paramVals(4,2,k,l,m)+...
                c_energy*paramVals(7,2,k,l,m);
            
            % distributions
            fX{1,k,l,m} = f_RD_Sa_DBE{k,l};
            fX{2,k,l,m} = f_RD_Sa_MCE{k,l};
            x0(1,k,l,m) = mean(f_RD_Sa_DBE{k,l});
            x0(2,k,l,m) = mean(f_RD_Sa_MCE{k,l});
            for i = 3:n_x
                if strcmp(pds{i,l},'lognormal')
                    % assumed params provided are median and dispersion
                    paramVals(i,1,k,l,m) = log(paramVals(i,1,k,l,m));
                end
                fX{i,k,l,m} = makedist(pds{i,l},paramNames{i,1,l},...
                    paramVals(i,1,k,l,m),...
                    paramNames{i,2,l},paramVals(i,2,k,l,m));
                x0(i,k,l,m) = mean(fX{i,k,l,m});
            end
            Rzz(:,:,k,l,m)   = eye(n_x);     % Start with no correlation
            Rzz(1,2,k,l,m)   = r_x1_x2(k,l); % Empirical; limited "IDA"
            Rzz(1,3,k,l,m)   = r_x1_x3(k,l); % Empirical; Monte Carlo
            Rzz(2,3,k,l,m)   = r_x2_x3(k,l); % Empirical; Monte Carlo
            Rzz(4,5,k,l,m)   = min([1,(1+ELMR+EL_PP_cases(k,l))*...
                paramVals(4,2,k,l,m)/paramVals(5,2,k,l,m)]);
            Rzz(7,5,k,l,m)   = c_energy*paramVals(7,2,k,l,m)/...
                paramVals(5,2,k,l,m);
            Rzz(:,:,k,l,m) = Rzz(:,:,k,l,m)+...
                transpose(Rzz(:,:,k,l,m))-eye(n_x);
            while any(eig(Rzz(:,:,k,l,m))<=0)
               disp(['Non-positive-definite correlation matrix for'...
                   ' struc, case, energy combo: ' num2str([l, k, m])]);
               [vec, val] = eig(Rzz(:,:,k,l,m));
               ind = vec(:,1)~=0;
               disp('Following variables cause problems & are adjusted:');
               disp(x(ind));
               eps = zeros(n_x);
               R_0 = Rzz(ind,ind,k,l,m)==0;
               R_neg = Rzz(ind,ind,k,l,m)<0;
               eps(ind,ind) = R_adj*(ones(sum(ind))-eye(sum(ind))-R_0);
               eps(ind,ind) = eps(ind,ind)-2*0.01*R_neg; % reduce negative correlation
               Rzz(:,:,k,l,m) =  Rzz(:,:,k,l,m) - eps;
               if(R_err>R_adj_max)||(max(max(abs(Rzz(:,:,k,l,m))))>1)
                  error('Unable to fix correlation by adjusting first eigenvector.')
                  break
               end
               if R_err==0
                  record_Rzz_err = [record_Rzz_err;k, l, m]; 
               end
               R_err = R_err + 1;
               R_err_max = max([R_err_max,R_err*R_adj]);
            end
            R_err = 0;
        end
    end
end
% Hard-code tough cases for HL-RF settings to be used in FORM.
x0(1:3,2,1,:) = [0.0112;0.0141;0.0439]*ones(1,numPPeVal);
x0(1:3,4,1,:) = [0.0176; 0.0238; 0.0908]*ones(1,numPPeVal); 
x0(1:3,1,2,:) = [0.0102;0.0177;0.043485]*ones(1,numPPeVal);
x0(1:3,6,2,:) = [0.0102;0.0177;0.043485]*ones(1,numPPeVal);
x0(1:3,2,2,:) = [0.0102;0.0141;0.0462]*ones(1,numPPeVal);
x0(1:3,3,2,:) = [0.0102;0.0141;0.047]*ones(1,numPPeVal);
lambda            = 0.15*ones(m_g,numCase,numStruct);
lambda(3,:,:)     = 0.05*ones(size(lambda(3,:,:)));
lambda(3,2,1)     = 0.001;
lambda(3,[3,6],2) = 0.001;
lambda(3,2,2)     = 0.0001;
lambda(3,1,2)     = 0.1;
%% 8.  Module 3: Decision support for individual limit states using FORM
pf        = zeros(m_g, numCase, numStruct, numPPeVal);
beta_FORM = zeros(m_g, numCase, numStruct, numPPeVal);
alpha     = zeros(n_x, m_g, numCase, numStruct, numPPeVal);
u_star    = zeros(n_x, m_g, numCase, numStruct, numPPeVal);
x_star    = zeros(n_x, m_g, numCase, numStruct, numPPeVal);
gamma     = zeros(n_x, m_g, numCase, numStruct, numPPeVal);

for m = 1:numPPeVal
    disp(['theta_envelope: ' PP_e(m),' ========'])
    for l = 1:numStruct
        disp(['Structure: ' strucs{l},' --------'])
        rd = mu_RD_c(l)*transpose(RD_c_fac);
        for k = 1:numCase
            disp(['   Case: ' num2str(k)])
            % set up g3: case-specific and nonlinear
            coefn = coeffnames(P_c_fit{k,l});
            cofv  = coeffvalues(P_c_fit{k,l});
            g3 = g3_formula;
            for i = 1:size(coefn,1)
                g3 = subs(g3,coefn{i},cofv(i));
            end
            grad_g3 = sym(zeros(n_x,1));
            grad_g3(3) = diff(g3,1);
            gs{3} = matlabFunction(g3,'vars',{x});
            grad_gs{3} = matlabFunction(grad_g3,'vars',{x});
            fX_i = reshape(fX(:,k,l,m),n_x,1);
            
            % loop over limit state functions; g3 special treatment
            for j = 1:m_g
                disp(['      g: ' num2str(j)])
                % can only use reliability approach if reasonable collapse
                % probabilities and some conditional probability greater 
                % than unacceptable value (0.1)
                if (j~=3)||((j==3)&&(P_c(k,l)>pC_min)&&...
                        (sum(P_c_fit{k,l}(rd)>=a(3,n_x+1))>0))
                    [pf(j,k,l,m), beta_FORM(j,k,l,m), alpha(:,j,k,l,m), ...
                        u_star(:,j,k,l,m), x_star(:,j,k,l,m), ...
                        gamma(:,j,k,l,m)] = M3_FORM(fX_i, ...
                        gs{j}, grad_gs{j},...
                        Rzz(:,:,k,l,m), x0(:,k,l,m),lambda(j,k,l),0);
                else
                    disp(['        Setting pf3 to 0'...
                                '  for [struc case PPe]: '...
                                num2str([l, k, m])...
                                ]);
                    record_g3_err = [record_g3_err;k, l, m];
                    pf(j,k,l,m) = 0;
                    beta_FORM(j,k,l,m) = 6;
                    alpha(:,j,k,l,m) = -inv(transpose(...
                        chol(Rzz(:,:,k,l,m))))*Rzz(:,3,k,l,m);
                    u_star(2:3,j,k,l,m) = [6;6];
                end
            end
        end
    end
end
%% 9.  Module 3: Systems-FORM reliability for cut sets
disp('Analyzing cut sets');
cuts{1} = 1;
cuts{2} = 3;
cuts{3} = 4;
cuts{4} =[5,6,7,8];
n_c = size(cuts,2);
pF_cut = zeros(n_c, numCase,numStruct,numPPeVal);
for m = 1:numPPeVal
    for l = 1:numStruct
        for k=1:numCase
            for I=1:n_c
                n_g = size(cuts{I},2);
                R_UU = eye(n_g);
                Beta = zeros(n_g,1);
                for ii = 1:(n_g-1)
                    i_c = cuts{I}(ii);
                    for jj = (ii+1):n_g
                        j_c = cuts{I}(jj);
                        alpha_i = alpha(:,i_c,k,l,m);
                        alpha_j = alpha(:,j_c,k,l,m);
                        R_UU(ii,jj) = transpose(alpha_i)*alpha_j;
                        R_UU(jj,ii) = R_UU(ii,jj);
                    end
                    Beta(ii) = beta_FORM(i_c,k,l,m);
                end
                Beta(n_g) = beta_FORM(cuts{I}(n_g),k,l,m);
                pF_cut(I,k,l,m) = mvncdf(-Beta,zeros(n_g,1),R_UU);
            end
        end
    end
end
%% 10. Module 3: Inclusion-Exclusion Rule for Final system reliability
% Probabilities of intersections of events, use upper + lower bounds
pF_upper1 = reshape(sum(pF_cut,1),numCase,numStruct,numPPeVal);

% intersection order = 2
disp('Computing intersection of 2 cut sets');
IJ = combnk(1:4,2);
nComb = size(IJ,1);
pF_cut2 = zeros(nComb,numCase,numStruct,numPPeVal);
pF_lower2 = zeros(size(pF_upper1));
for m = 1:numPPeVal
    for l = 1:numStruct
        for k=1:numCase
            for II = 1:nComb
                I = IJ(II,1);
                J = IJ(II,2);
                c = unique([cuts{I},cuts{J}]);
                n_g = size(c,2);
                R_UU = eye(n_g);
                Beta = zeros(n_g,1);
                for ii = 1:(n_g-1)
                    i_c = c(ii);
                    for jj = (ii+1):n_g
                        j_c = c(jj);
                        alpha_i = alpha(:,i_c,k,l,m);
                        alpha_j = alpha(:,j_c,k,l,m);
                        R_UU(ii,jj) = transpose(alpha_i)*alpha_j;
                        R_UU(jj,ii) = R_UU(ii,jj);
                    end
                    Beta(ii) = beta_FORM(i_c,k,l,m);
                end
                Beta(n_g) = beta_FORM(c(n_g),k,l,m);
                pF_cut2(II,k,l,m) = mvncdf(-Beta,zeros(n_g,1),R_UU);
            end 
            pF_lower2(k,l,m) = pF_upper1(k,l,m) - sum(pF_cut2(:,k,l,m));
        end
    end
end

% intersection order = 3
disp('Computing intersection of 3 cut sets');
IJK = combnk(1:4,3);
nComb = size(IJK,1);
pF_cut3 = zeros(nComb,numCase,numStruct,numPPeVal);
pF_upper3 = zeros(size(pF_upper1));
for m = 1:numPPeVal
    for l = 1:numStruct
        for k=1:numCase
            for III=1:nComb
                I = IJK(III,1);
                J = IJK(III,2);
                K = IJK(III,3);
                c = unique([cuts{I},cuts{J},cuts{K}]);
                n_g = size(c,2);
                R_UU = eye(n_g);
                Beta = zeros(n_g,1);
                for ii = 1:(n_g-1)
                    i_c = c(ii);
                    for jj = (ii+1):n_g
                        j_c = c(jj);
                        alpha_i = alpha(:,i_c,k,l,m);
                        alpha_j = alpha(:,j_c,k,l,m);
                        R_UU(ii,jj) = transpose(alpha_i)*alpha_j;
                        R_UU(jj,ii) = R_UU(ii,jj);
                    end
                    Beta(ii) = beta_FORM(i_c,k,l,m);
                end
                Beta(n_g) = beta_FORM(c(n_g),k,l,m);
                while any(eig(R_UU)<=0.0001)
                    disp(['Non-positive-definite limit state correlation matrix for'...
                        ' struc, case, energy combo: ' num2str([l, k, m])]);
                    [vec, val] = eig(R_UU);
                    ind = vec(:,1)~=0;
                    disp('Following limit states cause problems & are adjusted:');
                    disp(c(ind));
                    eps = zeros(n_g);
                    R_0 = R_UU(ind,ind)==0;
                    R_neg = R_UU(ind,ind)<0;
                    eps(ind,ind) = R_adj*(ones(sum(ind))-eye(sum(ind))-R_0);
                    eps(ind,ind) = eps(ind,ind)-2*R_adj*R_neg; % reduce negative correlation
                    R_UU =  R_UU - eps;
                    if(R_err>R_adj_max)||(max(max(abs(Rzz(:,:,k,l,m))))>1)
                        disp('Unable to fix correlation by adjusting first eigenvector.')
                        break
                    end
                    if R_err==0
                        record_R_UU_4_err = [record_R_UU_4_err;k, l];
                    end
                    R_err = R_err + 1;
                end
                R_err = 0;
                pF_cut3(III,k,l,m) = mvncdf(-Beta,...
                    zeros(n_g,1),R_UU);
            end
            pF_upper3(k,l,m) = pF_lower2(k,l,m) + sum(pF_cut3(:,k,l,m));
        end
    end
end

% intersection order = 4
disp('Computing intersection of 4 cut sets');
c = unique([cuts{1},cuts{2},cuts{3},cuts{4}]);
g_rank = [2;3;6]; % g3,4,7
pF_cut4 = zeros(numCase,numStruct,numPPeVal);
pF_lower4 = zeros(size(pF_upper1));
for m = 1:numPPeVal
    for l = 1:numStruct
        for k=1:numCase
            n_g = size(c,2);
            R_UU = eye(n_g);
            Beta = zeros(n_g,1);
            for ii = 1:(n_g-1)
                i_c = c(ii);
                for jj = (ii+1):n_g
                    j_c = c(jj);
                    alpha_i = alpha(:,i_c,k,l,m);
                    alpha_j = alpha(:,j_c,k,l,m);
                    R_UU(ii,jj) = transpose(alpha_i)*alpha_j;
                    R_UU(jj,ii) = R_UU(ii,jj);
                end
                Beta(ii) = beta_FORM(i_c,k,l,m);
            end
            Beta(n_g) = beta_FORM(c(n_g),k,l,m);
            while any(eig(R_UU)<=0.0001)
               disp(['Non-positive-definite limit state correlation matrix for'...
                   ' struc, case, energy combo: ' num2str([l, k, m])]);
               [vec, val] = eig(R_UU);
               ind = vec(:,1)~=0;
               disp('Following limit states cause problems & are adjusted:');
               disp(c(ind));
               eps = zeros(n_g);
               R_0 = R_UU(ind,ind)==0;
               R_neg = R_UU(ind,ind)<0;
               eps(ind,ind) = R_adj*(ones(sum(ind))-eye(sum(ind))-R_0);
               eps(ind,ind) = eps(ind,ind)-2*R_adj*R_neg; % reduce negative correlation
               R_UU =  R_UU - eps;
               if(R_err>R_adj_max)||(max(max(abs(Rzz(:,:,k,l,m))))>1)
                   disp('Unable to fix correlation by adjusting first eigenvector.')
                   break
               end
               if R_err==0
                   record_R_UU_3_err = [record_R_UU_3_err;k, l];
               end
               R_err = R_err + 1;
            end
            R_err = 0;
            pF_cut4(k,l,m) = mvncdf(-Beta,zeros(n_g,1),R_UU);
            pF_lower4(k,l,m) = pF_upper3(k,l,m) - pF_cut4(k,l,m);
        end
    end
end
pF = pF_lower4;
%% 11. Module 3: Ranking 
clc
[pF_rank, ord] = sort(pF(:));
alt_rank = names(ord);
disp(table(alt_rank(1:20,1),num2str(pF_rank(1:20)),...
    'VariableNames',{'Alt_config','P_F'}))
for m = 1:numPPeVal
    for l=1:numStruct
        pF_PPs(:,:,l,m) = [pF(1,l,m), pF(1,l,m), pF(1,l,m), pF(1,l,m);
            pF(2,l,m), pF(4,l,m), pF(6,l,m), pF(8,l,m);
            pF(3,l,m), pF(5,l,m), pF(7,l,m), pF(9,l,m)];
    end
end
%% 12. Module 2/3: Error diagnostics and save
%PF variance with inclusion/exclusion
save(['Data/',datestr(date,'yyyymmdd'),'_Flint_2019_M2_M3_all.mat']);
m = 3; % high-performing energy
for l=1:numStruct
    figure;
    hold on
    for k=1:numCase
        ps = [pF_upper1(k,l,m), pF_lower2(k,l,m), pF_upper3(k,l,m), pF_lower4(k,l,m)];
        plot(1:4,ps)
    end
end
%% 13. Module 3: Basic sensitivity
% first plotting example, rank 1,
l = 2; % steel
k = 1; % base case
m = 3; % high-perf
R_UU = eye(m_g);
for i = 1:(m_g-1)
    alpha_i = alpha(:,i,k,l,m);
    gammas(:,i) = gamma(:,i,k,l,m);
    xstars(:,i) = x_star(:,i,k,l,m);
    for j = (i+1):m_g
        alpha_j = alpha(:,j,k,l,m);
        R_UU(i,j) = transpose(alpha_i)*alpha_j;
        R_UU(j,i) = R_UU(i,j);
    end
end
for i=1:n_x
   mxs(i) =  mean(fX{i,k,l,m});
end
name1 = "steel, Base, PP_e: 1";
[mx wh] = max(strcmp(alt_rank,name1));
rank1 = wh;
pF_ex1 = pF(k,l,m);
Beta_ex1 = beta_FORM(:,k,l,m);
R_UU_ex1 = R_UU;
corr_j1 = 4;
pf_ex1 = pf(:,k,l,m);
gamma1 = gamma(:,corr_j1,k,l,m);
x_star1 = x_star(:,corr_j1,k,l,m);
pF_cut_ex1 = [pF_cut(1,k,l,m); NaN; pF_cut(2:4,k,l,m); NaN*ones(5,1)];
%disp(table(Beta_ex1,R_UU_ex1(:,corr_j)));
disp([name1,': rank: ',rank1,'; pF = ',pF_ex1])
disp(['Beta  |  pF_cut1  |  R_UU(:,' num2str(corr_j1),')']);
[Beta_ex1,pF_cut_ex1,R_UU_ex1(:,corr_j1)]

% second example, rank 2, conc d_pc = 0.5
l = 1; % conc
k = 6; % theta_pc 0.5
m = 3; % high-performing energy
R_UU = eye(m_g);
for i = 1:(m_g-1)
    for j = (i+1):m_g
        alpha_i = alpha(:,i,k,l,m);
        alpha_j = alpha(:,j,k,l,m);
        R_UU(i,j) = transpose(alpha_i)*alpha_j;
        R_UU(j,i) = R_UU(i,j);
    end
end
name2 = "conc, theta_pc: 0.5, PP_e: 1";
disp(name2)
[mx wh] = max(strcmp(alt_rank,name2));
rank2 = wh;
pF_ex2 = pF(k,l,m);
Beta_ex2 = beta_FORM(:,k,l,m);
R_UU_ex2 = R_UU;
corr_j2 = 6;
pf_ex2 = pf(:,k,l,m);
gamma2 = gamma(:,corr_j2,k,l,m);
x_star2 = x_star(:,corr_j2,k,l,m);
pF_cut_ex2 = [pF_cut(1,k,l,m); NaN; pF_cut(2:4,k,l,m); NaN*ones(5,1)];
%disp(table(Beta_ex1,R_UU_ex1(:,corr_j)));
disp([name2,': rank: ',rank2,'; pF = ',pF_ex2])
disp(['Beta  |  pF_cut1  |  R_UU(:,' num2str(corr_j2),')']);
[Beta_ex2,pF_cut_ex2,R_UU_ex2(:,corr_j2)]

% third example, concrete base, rank 26
l = 1; % conc
k = 1; % base
m = 3; % high-perf
R_UU = eye(m_g);
for i = 1:(m_g-1)
    for j = (i+1):m_g
        alpha_i = alpha(:,i,k,l,m);
        alpha_j = alpha(:,j,k,l,m);
        R_UU(i,j) = transpose(alpha_i)*alpha_j;
        R_UU(j,i) = R_UU(i,j);
    end
end
name3 = "conc, Base, PP_e: 1";
disp(name3)
[mx wh] = max(strcmp(alt_rank,name3));
rank3 = wh;
pF_ex3 = pF(k,l,m)
Beta_ex3 = beta_FORM(:,k,l,m);
R_UU_ex3 = R_UU;
corr_j3 = 3;
pf_ex3 = pf(:,k,l,m);
gamma3 = gamma(:,corr_j3,k,l,m);
x_star3 = x_star(:,corr_j3,k,l,m);
pF_cut_ex3 = [pF_cut(1,k,l,m); NaN; pF_cut(2:4,k,l,m); NaN*ones(5,1)];
%disp(table(Beta_ex1,R_UU_ex1(:,corr_j3)));
disp([name3,': rank: ',rank3,'; pF = ',pF_ex3])
disp(['Beta  |  pF_cut1  |  R_UU(:,' num2str(corr_j3),')']);
[Beta_ex3,pF_cut_ex3,R_UU_ex3(:,corr_j3)]
%% 14. Module 3: Theta/refined sensitivity
% select theta_pc, concrete and g5 (LCC EE)  for calculation
j = 5;
k = 6;
theta_k = theta(2);
k_next = 7; %dpc=1
k_base = 1;
ks = [k_base, k, k_next];
l = 1;
m = 3;
alpha_ex2 = alpha(:,j,k,l,m);
u_star2 = u_star(:,j,k,l,m);
x_star2 = x_star(:,j,k,l,m);
Rzz_ex2 = Rzz(:,:,j,l,m);
L = chol(Rzz_ex2)';
Linv = inv(L);
fZ = makedist('normal',0,1);
for kk=1:3
    z_star(:,kk) = L*u_star(:,j,ks(kk),l,m);
end
J_UZ = Linv;
J_Ztheta = zeros(n_x);
dXdTheta = zeros(n_x,1);
dZdTheta = zeros(n_x,1);
for i=1:n_x 
    dXdTheta(i) = mean(reshape(x_star(i,j,ks(2:3),l,m) -...
    x_star(i,j,ks(1:2),l,m),2,1)/0.5);
    dZdTheta(i) = mean((z_star(i,2:3) - z_star(i,1:2))/0.5);
    if any(i-[4,6,7]==0) % normally distributed with given CV
        J_Ztheta(i,i) = -x_star2(i)/(paramVals(i,1,k,l,m)^2*CV(i))*...
            pdf(fX{i,k,l,m},x_star2(i)) /  pdf(fZ,z_star(i));
    else %if i==5
            
        %else% use numerical derivative
            J_Ztheta(i,i) = 1;%dZdTheta(i);
        %end
    end
end
dMudTheta = [NaN; NaN; NaN; % obtaining from dZdtheta directly
    c_case(ks(3),l); % x4 slope = cost increase for PP
    NaN; % obtaining from dZdtheta directly
    e_case(ks(3),l);
    0];
dZdTheta([4,6:7]) = J_Ztheta([4,6:7],[4,6:7])*dMudTheta([4,6:7]);
dBdTheta = mean(reshape(beta_FORM(j,ks(2:3),l,m) -...
    beta_FORM(j,ks(1:2),l,m),2,1)/0.5);
dBetadTheta = alpha_ex2'*J_UZ*dZdTheta;
[dBetadTheta, dBdTheta]
%% 15. Write tabular output
% COST AND ENERGY
% concrete, steel
% lateral PP - cost_slope ee_slope 
m_l_out = dataset(PPnames(1:3)', c_PP_l(:,2,1), e_PP_l(:,2,1),...
    c_PP_l(:,2,2), e_PP_l(:,2,2),'VarNames',{'PP','a_cost_c','a_ee_c',...
    'a_cost_s','a_ee_s'});
export(m_l_out,'file','Data/m_lateral_slopes.txt','Delimiter','tab');

% envelope PP cost_slope ee_slope
m_e_out = dataset({'window/wall ratio'}, c_PP_e(3,1) - c_PP_e(2,1),...
  e_PP_e(3,1) - e_PP_e(2,1),'VarNames',{'PP','a_cost','a_ee'});
export(m_e_out,'file','Data/m_envelope_slopes.txt','Delimiter','tab');

% economic and environmental ranges
% concrete - x4 x5 x6 x7 central value range across k,m
% steel - x4 x5 x6 x7 central value range across k,m
ranges = cell(2,4);
dash = ['-';'-';'-';'-'];
for l=1:numStruct
    vals = reshape(paramVals(4:7,1,:,l,:),4,numCase*numPPeVal);
    for i=1:4
       ranges{l,i} = [num2str(round(min(vals(i,:),[],2),3,'significant')),...
       '-',...
       num2str(round(max(vals(i,:),[],2),3,'significant'))];
    end
   
end
ranges
[mins whmin] = min(vals(3,:),[],2);
[maxs whmax] = max(vals(3,:),[],2);
namesNS{[whmin, whmax]}

% STRUCTURAL PERFORMANCE
% concrete, steel
% config - a b beta a b beta
PBEE_param_out = dataset(caseNames,paramNCs(1,:,1)', paramNCs(2,:,1)',...
    paramNCs(3,:,1)', paramNCs(1,:,2)', paramNCs(2,:,2)',...
    paramNCs(3,:,2)', 'VarNames',{'Config','a conc','b conc',...
    'beta conc','a steel','b steel','beta steel'});
export(PBEE_param_out,'file','Data/PBE_median_drift.txt','Delimiter','tab');

% concrete, steel
% config medX1 betaX1 medX2 betaX2 medX3 sigmaX3 rhoX1X2 rhoX2X3...
for l=1:numStruct
    for k=1:numCase
        medX1(k,l) = f_RD_Sa_DBE{k,l}.mu;
        betaX1(k,l) = f_RD_Sa_DBE{k,l}.sigma;
        medX2(k,l) = f_RD_Sa_MCE{k,l}.mu;
        betaX2(k,l) = f_RD_Sa_MCE{k,l}.sigma;
        medX3(k,l) = pd_Cs{k,l}.mu;
        betaX3(k,l) = pd_Cs{k,l}.sigma;
        
    end
end
X1X2X3_param_out = dataset(caseNames, medX1(:,1), betaX1(:,1), medX2(:,1),...
    betaX2(:,1), medX3(:,1), betaX3(:,1), r_x1_x2(:,1), r_x2_x3(:,1),...
    medX1(:,2), betaX1(:,2), medX2(:,2),...
    betaX2(:,2), medX3(:,2), betaX3(:,2), r_x1_x2(:,2), r_x2_x3(:,2),...
    'VarNames',{'Config','medX1 c','betaX1 c', 'medX2 c','betaX2 c', ...
        'medX3 c','betaX3 c', 'rhoX1X2 c','rhoX2X3 c',...
        'medX1 s','betaX1 s', 'medX2 s','betaX2 s', ...
        'medX3 s','betaX3 s', 'rhoX1X2 s','rhoX2X3 s'});
export(X1X2X3_param_out,'file','Data/X1X2X3_param.txt','Delimiter','tab');

% concrete, steel
% config pC_conv pf3 pC_conv pf3
collapse_out = dataset(caseNames, P_c(:,1), pf(3,:,1,1)',P_c(:,2),...
    pf(3,:,2,1)','VarNames',{'Config','Pc_conv_c','pf3_c','Pc_conv_s','pf3_s'});
export(collapse_out,'file','Data/collapse.txt','Delimiter','tab');

% ELLR AND OPERATIONAL PERFORMANCE
% concrete, steel
% ELLR median ELLR convolved
ellr_out = dataset(caseNames, EL_PP_cases(:,1),ELLR_c(:,1),...
    EL_PP_cases(:,2),ELLR_c(:,2),'VarNames',{'Config','ELLR_med_c',...
    'ELLR_conv_c','ELLR_med_s','ELLR_conv_s'});
export(ellr_out,'file','Data/ellr.txt','Delimiter','tab');

% concrete, steel
op_out = dataset(PP_e([2,1,3])', e_op(1,[2,1,3])',e_op(2,[2,1,3])',...
    'VarNames',{'PP_e','e_op_c','e_op_s'});
export(op_out,'file','Data/op_energy.txt','Delimiter','tab');

% FINAL RANKING/RESULTS
% write rank
rank_out = dataset(alt_rank,pF_rank);
export(rank_out,'file','Data/rank.txt','Delimiter','tab');

% write comparison of 3 examples
ex_out_pf = dataset([name1 name2 name3; rank1 rank2 rank3;...
    pF_ex1 pF_ex2 pF_ex3; corr_j1 corr_j2 corr_j3]);
export(ex_out_pf,'file','Data/examples_pF.txt','Delimiter','tab');
ex_out = dataset([Beta_ex1,pF_cut_ex1,R_UU_ex1(:,corr_j1),...
    Beta_ex2,pF_cut_ex2,R_UU_ex2(:,corr_j2),...
    Beta_ex3,pF_cut_ex3,R_UU_ex3(:,corr_j3)]);
export(ex_out,'file','Data/examples.txt','Delimiter','tab');

% write sensitivity of top-ranked example
js = [1,4:8];
xstars = xstars(:,js);
gammas = gammas(:,js);
ex_out_highest = dataset(mxs',...
    xstars(:,1), ...
    xstars(:,2), ...
    xstars(:,3), gammas(:,3), ...
    xstars(:,4), gammas(:,4), ...
    xstars(:,5), gammas(:,5), ...
    xstars(:,6), gammas(:,6));
export(ex_out_highest,'file','Data/examples_highest.txt','Delimiter','tab');
%% 16. Plot backbone curves
M2_Plot_Backbones
%% 17. Scatterplots of EDP versus IM
xlim = [0.015, 1.7];
ylim = [0.0005,0.6];
% xlim = [-4.5, 0.5];
% ylim = [-8.5,-0.5];
close all
for l=1:numStruct
    struc = strucs{l};
    RD_c  = mu_RD_c(l);
    xlab  = xlabs(l);
    
    % the ones for Fy
    M2_Plot_EDP_IM(xlim,ylim,Sas(:,l),EDPs(:,[1,4,5],l),RD_c,color2,...
        'F_y',struc,true,xlab,SAVE_FIG)
    % the ones for thetaP
    M2_Plot_EDP_IM(xlim,ylim,Sas(:,l),EDPs(:,[1,6,7],l),RD_c,color3,...
        'theta_p',struc,true,xlab,SAVE_FIG)
    %  the ones for thetapc
    M2_Plot_EDP_IM(xlim,ylim,Sas(:,l),EDPs(:,1:3,l),RD_c,color1,...
        'theta_pc',struc,true,xlab,SAVE_FIG)
    % the ones for Combined
    M2_Plot_EDP_IM(xlim,ylim,Sas(:,l),EDPs(:,[1,8,9],l),RD_c,color4,...
        'combined',struc,true,xlab,SAVE_FIG)

end
%% 18. Plot fragility curves for damage states
xlims = [0,0.2; 0, 0.5; 0, 1; 0, 2; 0, 2];
for l=1:numStruct
    struc = strucs{l};
    RD_c  = mu_RD_c(l);
    xlab  = xlabs(l);
    Lstyle = Lstyles{l};
    for k = 1:numPP % PP types
        for kk = 1:numDS
            M2_Plot_Fragility(xlims(kk,:),SaCalc,...
                reshape(fragSs(:,kk,PPcols(l,:),l),length(SaCalc),3),...
                Lstyle, colors(k,:), DSnames{kk},...
                [PPnames{k}, num2str(RD_c)],struc,true,true,true)
        end
        M2_Plot_Fragility(xlims(numDS+1,:),SaCalc,...
            reshape(DamageProbs(:,numDS+1,PPcols(l,:),l),length(SaCalc),3),...
            Lstyle, colors(k,:), DSnames{5},...
            [PPnames{k},num2str(RD_c)],struc,true,true,true)
    end
end

for j = 1:numPP % PP types
    for kk = 1:(numDS + 1)
        name = ['frag_',DSnames{kk},'_',PPnames{j},...
            num2str(mu_RD_c(1)),'_conc.fig'];
        figconc = openfig(name);
        name = ['frag_',DSnames{kk},'_',PPnames{j},...
            num2str(mu_RD_c(2)),'_steel.fig'];
        figsteel = openfig(name);
        L = findobj(figconc,'type','line');
        copyobj(L,findobj(figsteel,'type','axes'));
        if SAVE_FIG==true
            savefig([pwd,'/Figs/fig/frag_',DSnames{kk},'_',PPnames{j},...
                '_RD_c_1p5_median_conc_steel.fig']);
        end
        print('-depsc', [pwd,'/Figs/eps/frag_',DSnames{kk},...
            '_',PPnames{j},'_RD_c_1p5_median_conc_steel.eps']);
    end
end
%% 19. Scatterplot of ELLR values; concrete and steel together
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
    'YMinorTick','on','YTickLabel',{'0.1%','0.3%','0.5%','0.7%','0.9%'},...
    'units','inches','position',[0.57 0.42 2.35 2.45]);
box(gca,'on');
hold(gca,'all');
MstylesComb = {'o';'s';'^';'V'};
for l = 1:numStruct
    Lstyle = Lstyles{l};
    for k = 1:numPP
        Mstyle = MstylesComb{k};
        plot(theta,EL_PPs_c(:,k,l),Lstyle,'MarkerEdgeColor',colors(k,:),...
            'Marker',Mstyle,...
             'Color',colors(k,:))
    end
end
xlabel('Performance Parameter');
ylabel('Expected Lifetime Loss Ratio');
if SAVE_FIG==true
    savefig([pwd,'/Figs/fig/EL_conc_steel_RD_c_all.fig']);
end
print('-depsc', [pwd,'/Figs/eps/EL_conc_steel_RD_c_all.eps']);
%% 20. Plot of final results for each configuration
gcf = figure('Color',[1 1 1]);
set(gcf, 'units','inches','position',[1 1 3 3],'PaperUnits', 'Inches');
gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
            'FontSize',10,...
    'FontName','Arial',...
    'Linewidth', 1,...
    'TickLength', [0.02 0.035],... 
    'YLim',[0, 1],'YMinorTick','on',...
    'XTick',[0,0.5,1],'Xcolor',[0,0,0],'Ycolor',[0,0,0],...
    'units','inches','position',[0.57 0.42 2.35 2.45]);
box(gca,'on');
hold(gca,'all');
for m = [3,2,1]
    Lwidth = Lwidths(m);
    for l = 1:numStruct
        Lstyle = Lstyles{l};
        for k = 1:numPP
            Mstyle = MstylesComb{k};
            plot(theta,pF_PPs(:,k,l,m),Lstyle,'LineWidth',Lwidth,...
            'MarkerEdgeColor',colors(k,:), 'Marker',Mstyle,...
                'Color',colors(k,:))
        end
    end
end
xlabel('Performance Parameter Value');
ylabel('Probability of System Failure');
if SAVE_FIG==true
    savefig([pwd,'/Figs/fig/pF_all.fig']);
end
print('-depsc', [pwd,'/Figs/eps/pF_all.eps']);
%% 21. Pareto-ish front of 3 cut set failure probabilities
 gcf = figure('Color',[1 1 1],'units','inches','position',... %'none'
     [1 1 4 4],'PaperUnits', 'Inches');
   gca = axes('Parent',gcf,'YGrid','off','XGrid','off',...
            'FontSize',10,...
    'FontName','Arial',...
    'Linewidth', 1,...
    'TickLength', [0.02 0.035],... 
    'XLim',[0,1],'XTick',[0,0.2,0.4,0.6,0.8,1],'XColor','black',...
    'YLim',[0,1],'YTick',[0,0.2,0.4,0.6,0.8,1],'YColor','black',...
    'ZLim',[0,1],'ZTick',[0,0.2,0.4,0.6,0.8,1],'ZColor','black');
 box(gca,'on');
 grid on;
 hold(gca,'all');
 daspect([1 1 1]);
 e_cols = [0.5, 0.5, 0.5;
           0.75, 0.75, 0.75;
           0.25,0.25,0.25];
 for m = 1:numPPeVal
     for l = 1:numStruct
         x = reshape(pF_cut(2,:,l,m),numCase,1);
         y = reshape(pF_cut(3,:,l,m),numCase,1);
         z = reshape(pF_cut(4,:,l,m),numCase,1);
         e_col = e_cols(m,:);
         for k = 1:4
             Mstyle = MstylesComb{k};
             if l==1
                 plot3(x(1),y(1),z(1),...
                     'MarkerFaceColor',colors(k,:),...
                     'Marker',Mstyle,...
                     'MarkerSize',4,...
                     'MarkerEdgeColor',e_col)
                 plot3(x(2*k),y(2*k),z(2*k),...
                     'MarkerFaceColor',colors(k,:),...
                     'Marker',Mstyle,...
                     'MarkerSize',8,...
                     'MarkerEdgeColor',e_col)
                 plot3(x(2*k+1),y(2*k+1),z(2*k+1),...
                     'MarkerFaceColor',colors(k,:),...
                     'Marker',Mstyle,...
                     'MarkerSize',12,...
                     'MarkerEdgeColor',e_col)
             else
                 plot3(x(1),y(1),z(1),...
                     'MarkerFaceColor',e_col,....
                     'Marker',Mstyle,...
                     'MarkerSize',4,...
                     'MarkerEdgeColor',colors(k,:))
                 plot3(x(2*k),y(2*k),z(2*k),...
                     'MarkerFaceColor',e_col,...
                     'Marker',Mstyle,...
                     'MarkerSize',8,...
                     'MarkerEdgeColor',colors(k,:))
                 plot3(x(2*k+1),y(2*k+1),z(2*k+1),...
                     'MarkerFaceColor',e_col,...
                     'Marker',Mstyle,...
                     'MarkerSize',12,...
                     'MarkerEdgeColor',colors(k,:))
             end
         end
     end
 end
%view([-37.5 30]);
view([-77.5 16.3]);
xlabel('Pf II: Collapse');
ylabel('Pf III: Initial Cost');
zlabel('Pf IV: Tradeoffs');
if SAVE_FIG==true
    savefig([pwd,'/Figs/fig/pF_II_III_IV.fig']);
end
print('-depsc', [pwd,'/Figs/eps/pF_II_III_IV.eps']);
%% 22. Plots of 3D contours and limit state functions: not used in paper
close all
clc
m = 3;
k = 1;
l = 2;
coefn = coeffnames(P_c_fit{k,l});
cofv  = coeffvalues(P_c_fit{k,l});
g3 = g3_formula;
for i = 1:size(coefn,1)
    g3 = subs(g3,coefn{i},cofv(i));
end
gs{3} = matlabFunction(g3,'vars',{x});

% last cut set = 4 limit state functions
is = [1,2,3]; % random variables to plot
js = [1,2,3,9,10]; % limit state functions to plot
M3_Plot_Hyperpolygon(gs,fX(:,k,l,m),Rzz(:,:,k,l,m),is,js,...
    u_star(:,:,k,l,m)',alpha(:,:,k,l,m)',beta_FORM(:,k,l,m),...
    configNames{k,l},1,1,0,0);
%% End