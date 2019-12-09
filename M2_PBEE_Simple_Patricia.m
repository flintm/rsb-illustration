%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code imports EDP and IM Results for two different SDOF systems and
% calculate their expected loss
% Mohsen Zaker Esteghamati, modified by Madeleine Flint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sa, EDP, pd_C, paramNC, fragS, DamageProb, SaCalc, EL, EAL_S, EAL_NS, EAL_C] = M2_PBEE_Simple_Patricia(struc, EDP_c, Years)
%%%Note%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Throughout this code the data are arranged as matrices by 9 columns to
% represent 9 cases under study:
%   First column: original structure values
%   Columns 2,3: theta_pc multiplied by 1.5 and 2
%   Columns 4,5: F_y multiplied by 1.5 and 2
%   Columns 6,7: theta_p multiplied by 1.5 and 2
%   Columns 8,9: Combined effect of all PPs multiplied by 1.5 and 2
% This organization does not match that used in the paper; the 
% M2_M3_Main_Run_All function handles the re-ordering.

%% Load data for structure under study
%number of cases under study
numCase=1;
%importing Opensees results saved as matlab variables
if strcmp(struc,'steelMod')
    %spectral acceleration (normalized by g)
    Sa = importdata('Data/Sa_steelMod.mat');
    %roof drift values
    EDP = importdata('Data/EDP_steelMod.mat');
else
    Sa = importdata('Data/Sa_masonry.mat');
    EDP = importdata('Data/EDP_masonry.mat');

end
%% Collapse and Non-Collapse Analysis
% getting collapse indices
col_ind = EDP >= EDP_c;

% zero matrix for regression coefficient for non-collapse data
paramNC = zeros(3,numCase);

% parameters of noncollapse model
for i=1:numCase
    % fit a regression to non-collapse data and sort them in a 3X9
    % matrix where each row represents a parameter and each column 
    % represents a case 
    md=fitlm(log(Sa(~col_ind(:,i))),log(EDP(~col_ind(:,i),i)));
    paramNC(1,i) = md.Coefficients{1,1}; % intercept
    paramNC(2,i) = md.Coefficients{2,1}; % slope
    paramNC(3,i) = md.RMSE;              % beta
end

% Fit a log-logistic regression to collapse data
% Sa interval to evaluate the logistic distribution
SaCalc = (0:0.001:2)';

%creating the zero array and cell
logitFit = zeros(length(SaCalc),numCase);
pd_C     = cell(numCase,1);
%now looping over each case and data
for i=1:numCase
    if sum(col_ind(:,i))>1
        pd_C{i,1} = fitdist(log(Sa(col_ind(:,i))),'logistic');
        logitFit(:,i)=cdf(pd_C{i,1},log(SaCalc)); 
    else
        if sum(col_ind(:,i))==1
            sa_c = Sa(col_ind(:,i));
            logitFit(SaCalc>=sa_c,i) = 1; % empirical CDF
            % but still need to make a distribution for collapse
            % uncertainty
            pd_C{i,1} = makedist('logistic','mu',log(sa_c),'sigma',0.01);
        else % no collapse recordings
            % logitfit already all 0s, but still need distribution
            pd_C{i,1} = makedist('logistic','mu',log(5),'sigma',0.01); 
        end
    end
end
%% Fragilities
% Limit states based on HAZUS high-code, using roof drifts
% for structural and nonstructural damage states
fragS=zeros(length(SaCalc),4,1); 
fragNS=zeros(length(SaCalc),4,1); 
% defining the limit states drifts
if strcmp('struc','conc') % **
        DSS=[1.5,3,9,24]/648; % concrete, 648 = roof height, inches
        DSNS=[1.80,3.60,11.25,22.50]/648; % concrete
else
        DSS=[2.16,4.32,10.80,28.80]/612;% steel, 612 = height, inches
        DSNS=[2.16, 4.32, 13.50, 27.0]/612; %steel
end

numDS = size(DSS,2);
numNSDS = size(DSNS,2);

% Fragility formulation for cloud analysis 
for i=1:numCase
    for j=1:numDS
        % j is the indice for damage state and i is the one for case study
        % (thetap*1.5,thepc*1.5,etc)
        fragS(:,j,i) = 1-normcdf((log(DSS(j))-paramNC(1,i)-...
            paramNC(2,i)*log(SaCalc))/paramNC(3,i));
        fragNS(:,j,i) = 1-normcdf((log(DSNS(j))-paramNC(1,i)-...
            paramNC(2,i)*log(SaCalc))/paramNC(3,i));
    end
end
%% Mean annual frequency of exceeding damage

% hazardinfo is taken from USGS where first column is Sa values (in g) and 
% the second and the third column shows exceedance of SA values for periods 
% of 1 and 2 seconds. 
load('Data/hazardInfo.mat')
    if strcmp('struc','conc')% **
        T_interp=1.16; % period of concrete structure **
    else
        T_interp=1.32; % period of concrete structure **
    end
% Linear interpolation in log-log space is used for building period.
 for i=1:length(hazardInfo)
         hazardInterp(i)=interp1([1,2],[log(hazardInfo(i,3)),...
             log(hazardInfo(i,4))],T_interp,'linear');
end

% Reformat as vector with SA and exceedance as usual; get dLambda.
hazardCurve=[hazardInfo(:,1),exp(hazardInterp)'];
h=0.0000001;
for k = 1:length(SaCalc)
    % Interpolating value of F(x+h)
    fxp = interp1(hazardCurve(:,1),hazardCurve(:,2),SaCalc(k)+h,'spline'); 
    %  Interpolating value of F(x-h)
    fxn = interp1(hazardCurve(:,1),hazardCurve(:,2),SaCalc(k)-h,'spline');  
    % Two-point formula to compute numerical differentiation
    PDF_temp(k) = abs((fxp - fxn))/(2*h);             
end                                      
dlambdaIM = smooth(PDF_temp);
%% Damage probabilities
%DamageProb is the probability that structure is within a damage state
DamageProb=zeros(length(SaCalc),numDS + 1 ,numCase);
DamageProbNS=zeros(length(SaCalc),numDS + 1 ,numCase);
for i=1:numCase
    for k=1:length(SaCalc)
        for j=2:numDS % (pre-collapse)
            DamageProb(k,j-1,i) = (fragS(k,j-1,i)-...
                fragS(k,j,i)).*(1-logitFit(k,i));
            DamageProbNS(k,j-1,i) = (fragNS(k,j-1,i)-...
                fragNS(k,j,i)).*(1-logitFit(k,i));
        end
        DamageProb(k,4,i) = fragS(k,4,i).*(1-logitFit(k,i));
        DamageProbNS(k,4,i) = fragNS(k,4,i).*(1-logitFit(k,i));
        DamageProb(k,5,i) = logitFit(k,i); % collapse for structural only
     end
end
%% HAZUS commercial buildings repair cost under structural damage states
% as % of total replacement cost
Repair_cost=[0.4,1.9,9.6,19.2,110]; % 4 DS + collapse *
Repair_cost_NS=[0.7 ,3.3, 16.4, 32.9]; % 4 DS
% Convolution to obtain Mean Annual Frequency of Exceeding repair cost
for i = 1:numCase
    for kk = 1:(numDS + 1)    
        for j=1:length(dlambdaIM)
             conv(j,kk,i) = dlambdaIM(j).*DamageProb(j,kk,i);
        end
        % this is the mean annual frequency of exceeding each damage state
        MAF(kk,i)=trapz(SaCalc,conv(:,kk,i));
    end
    for kk = 1:(numNSDS)    
        for j=1:length(dlambdaIM)
             conv_NS(j,kk,i) = dlambdaIM(j).*DamageProbNS(j,kk,i);
        end
        % this is the mean annual frequency of exceeding each damage state
        MAF_NS(kk,i)=trapz(SaCalc,conv_NS(:,kk,i));
    end 
end

% Multiplying MAFE by expected repair cost and obtaining lifetime
% distributions
for i=1:numCase
    for kk = 1:(numDS + 1) 
    EAL_DS(kk,i)=Repair_cost(kk)*MAF(kk,i);
    end
    for kk = 1:(numNSDS)    
       EAL_NS(kk,i)=Repair_cost_NS(kk)*MAF_NS(kk,i);
    end
end
EAL_S = sum(EAL_DS(1:4,:));
EAL_C = EAL_DS(5,:);
EAL_NS = sum(EAL_NS);
EAL = EAL_S + EAL_C + EAL_NS;

%  Multiplying EAL by structure life time to get the expected loss in Years of structure
EL = EAL*Years;
end

