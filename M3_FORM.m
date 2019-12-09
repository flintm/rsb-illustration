% FORM calculations using improved HLRF algorithm, assuming Nataf distrib
% with correlated RVs.
% For Flint et al. 2019, written by Madeleine Flint, 2019-09-19
function [pf, beta_FORM, alpha, u_star, x_star, gamma] = M3_FORM(fX, g,grad_g, Rzz, x0, lambda, FLAG_PLOT)
% fX     = cell array of probability distribution objects
% g      = function handle for limit state function g(X)
% grad_g = function handle for gradient of g(X)
% R_zz   = assumed correlation between Nataf variables in Z-space
%% Set up FORM tolerances/controls
it_max = 400;
eps_1 = 1e-4;
eps_2 = 1e-4;
lambda_min = 1e-06;
%% Nataf transform and initial guess setup
n = size(Rzz,1);
z0 = zeros(n,1);

fZ = makedist('normal',0,1);
L = chol(Rzz)';
Linv = inv(L);

for i=1:n
    z0(i) = icdf(fZ, cdf(fX{i},x0(i)));
end
x = x0;
z = z0;
z(z0<-8) = -8;
z(z0>8) = 8;
u = Linv*z; 
beta = norm(u);
J_ZX = zeros(n);
%% Iteration: improved-Hasofer-Lind-Rackwitz-Feissler
done=0;
k = 1;

while(~done)
    % compute the variables needed for one iteration
    g_x = g(x(:,k));
    h_u = g_x;
    
    % Jacobian of h(U)
    for i=1:n
        J_ZX(i,i) = pdf(fX{i},x(i,k)) / pdf(fZ,z(i,k));
    end
    J_UX = Linv*J_ZX;
    J_XU = inv(J_UX);
    
    grad_h = J_XU' * grad_g(x(:,k));
    norm_h = norm(grad_h);
    alpha = - grad_h/norm_h;
    beta_old = beta(k);
    beta(k+1) = alpha' * u(:,k);
    
    % calculate new estimates of u, z, and x
    u(:,k+1) = alpha * (beta(k+1) + lambda*h_u/norm_h);
    % too far out of bounds, reduce lambda
    while (sum(abs(u(:,k+1))>9)>0)||(sum(isnan(u(:,k+1)))>0)
        disp('Reducing lambda');
        lambda = lambda*0.1;
        if lambda < lambda_min
            warning('Ended by NaN or Inf u for g_j(X)= %0.10s', char(g));
            done = 1;
            break
        end
        u(:,k+1) = alpha * (beta(k+1) + lambda*h_u/norm_h);
    end
    z(:,k+1) = L * u(:,k+1);
    for i=1:n
        x(i,k+1) = icdf(fX{i}, cdf(fZ,z(i,k+1)));
    end
    if (k >= 1) % check to see if we have converged
        check_1 = abs(h_u) < eps_1;
        check_2 = abs( beta(k+1) - beta_old) < eps_2;
        if (check_1 && check_2)
            done = 1; % we have converged to the design point
        end
    end
    if (k > it_max)
        warning('Ended by max iterations for g_j(X)= %0.10s', char(g))
        done = 1;
    end
    if (sum(isnan(u(:,k+1))|isinf(u(:,k+1)))>0)
       warning('Ended by NaN or Inf u for g_j(X)= %0.10s', char(g));
       done = 1;
    end
    k = k+1;
    %disp('**');
end

%% Record results and plot convergence (optional)
beta_FORM = beta(k);
x_star = x(:,k);
u_star = u(:,k);
pf     = normcdf(-beta_FORM);
Jcov   =  J_XU*J_XU';
Dp     = diag(sqrt(diag(Jcov))) ;
gamma  =  Dp*inv(J_XU)'*alpha;
gamma  = gamma/norm(gamma);

if(FLAG_PLOT)  
    figure
    plot(1:k, beta);
    xlabel('Iteration Number');
    ylabel('Reliability Index, beta');
    title(char(g))
end
end