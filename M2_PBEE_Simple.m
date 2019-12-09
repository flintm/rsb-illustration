%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code imports EDP and IM Results for two different SDOF systems and
% calculate their expected loss
% Mohsen Zaker Esteghamati, modified by Madeleine Flint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sa, EDP, pd_C, paramNC, fragS, DamageProb, ...
    EL_PP, EL_case, EAL_S, EAL_NS, EAL_C] = M2_PBEE_Simple(numCase,...
    EDP_c, hazardInfo, SaCalc, EDP, Sa, DSS, DSNS, Repair_cost_S, ...
    Repair_cost_NS, Repair_cost_C, T1,  Years)
%%%Note%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the illustration the data are arranged as matrices by 9 columns to
% represent 9 cases under study:
%   First column: original structure values
%   Columns 2,3: theta_pc multiplied by 1.5 and 2
%   Columns 4,5: F_y multiplied by 1.5 and 2
%   Columns 6,7: theta_p multiplied by 1.5 and 2
%   Columns 8,9: Combined effect of all PPs multiplied by 1.5 and 2
% This organization does not match that used in the paper; the 
% M2_M3_Main_Run_All function handles the re-ordering with the exception
% of EL_PP.

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
            warning('Only one collapse');
            pd_C{i,1} = makedist('logistic','mu',log(sa_c),'sigma',0.01);
        else % no collapse recordings
            % logitfit already all 0s, but still need distribution
            pd_C{i,1} = makedist('logistic','mu',log(5),'sigma',0.01); 
            warning('No collapses');
        end
    end
end
%% Fragilities
numDS = size(DSS,2);
numNSDS = size(DSNS,2);
if numDS~=numNSDS
   error('number of structural and nonstructural DS should be the same'); 
end
% for structural and nonstructural damage states
fragS=zeros(length(SaCalc),numDS,numCase); 
fragNS=zeros(length(SaCalc),numDS,numCase); 

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
%% Mean annual frequency of IM interpolation
% hazardinfo is taken from USGS where first column is Sa values (in g) and 
% the second and the third column shows exceedance of SA values for periods 
% of 1 and 2 seconds. 
% Linear interpolation in log-log space is used for building period.
 for i=1:length(hazardInfo)
         hazardInterp(i)=interp1([1,2],[log(hazardInfo(i,3)),...
             log(hazardInfo(i,4))],T1,'linear');
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
DamageProb   = zeros(length(SaCalc),numDS + 1, numCase);
DamageProbNS = zeros(length(SaCalc),numDS, numCase);
for i=1:numCase
    for k=1:length(SaCalc)
        for j=2:numDS % (pre-collapse)
            DamageProb(k,j-1,i) = (fragS(k,j-1,i)-...
                fragS(k,j,i)).*(1-logitFit(k,i));
            DamageProbNS(k,j-1,i) = (fragNS(k,j-1,i)-...
                fragNS(k,j,i)).*(1-logitFit(k,i));
        end
        DamageProb(k,numDS,i) = fragS(k,numDS,i).*(1-logitFit(k,i));
        DamageProbNS(k,numNSDS,i) = fragNS(k,numNSDS,i).*(1-logitFit(k,i));
        DamageProb(k,numDS+1,i) = logitFit(k,i); % collapse for structural only
     end
end
%% HAZUS commercial buildings repair cost under structural damage states
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
    for kk = 1:numDS 
    EAL_S(kk,i)=Repair_cost_S(kk)*MAF(kk,i);
    end
    for kk = 1:(numNSDS)    
       EAL_NS(kk,i)=Repair_cost_NS(kk)*MAF_NS(kk,i);
    end
    EAL_C(i) = Repair_cost_C*MAF(numDS+1,i);
end
EAL_S = sum(EAL_S);
EAL_NS = sum(EAL_NS);
EAL = EAL_S + EAL_NS + EAL_C;

%% reorganize to order desired for plotting (Fy, theta_p, theta_pc, all)
if numCase==9
    EAL_PP = [EAL(1), EAL(4:5); EAL(1), EAL(6:7); EAL(1), EAL(2:3); EAL(1), EAL(8:9);]';
    EAL    = EAL([1,4,5,6,7,2,3,8,9]);
else
    EAL_PP = EAL;
end

%  Multiplying EAL by structure life time Years
EL_PP = EAL_PP*Years;
EL_case = EAL*Years;
end

