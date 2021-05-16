%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2021 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
max_compliance = 10.0;
param.load = 10.0;
eeee = 2.0;
%
param.terminate = 10^(-7);
% param.num_data  = 1*10^(4);
% param.num_data2 = 1*10^(6);
param.num_data  = 1*10^(2);
param.num_data2 = 1*10^(2);
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data of the structure
% --->
[dll,matH,coord_x,ir,irr,ird] = member(0);
%
nk = size(coord_x,1); num.node   = nk;
nd = size(matH,1);    num.degree = nd;
nm = size(matH,2);    num.member = nm;
% <---
% Data of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
Idx.deg_of_load = 2;
%
vec_p = sparse(nd,1);
vec_p(Idx.deg_of_load) = -param.load;
% <---
% Load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix
% --->
sqrtK = matH * sparse( diag(sqrt(eeee ./ dll)) );
% <---
% Stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nominal optimization
% --->
fprintf(' =========================================================== \n');
fprintf(' ==  Nominal optimization >>> \n');
cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    variable nominal_x(nm,1)
    minimize( dll' * nominal_x )
    subject to
        nominal_x >= 0;
        [max_compliance, vec_p';
            vec_p, (sqrtK * diag(nominal_x) * sqrtK')] >= 0;
cvx_end

nom_matK = sqrtK * diag(nominal_x) * sqrtK';
fprintf('     compliance = %6.5f [J] ; bound = %6.5f [J} \n',...
    (vec_p' * (nom_matK \ vec_p) * 10) , (max_compliance * 10) );
fprintf('     opt.val. = %3.5e [cm^3] = %3.5e [mm^3] \n',...
    (dll' * nominal_x * 100), (dll' * nominal_x * (10^5)) );

fprintf(' =========================================================== \n');
% <---
% Nominal optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robust RBDO  optimization
% --->
param.alpha = 0.2;
param.beta  = 0.01;
param.eps = 0.01;
tilde_mu = 0 * ones(nm,1);
tilde_Sigma = (0.05 * eye(nm)) + (0.02 * ones(nm,nm));
param_kappa    = -norminv(param.eps);
param_kappa_sq =  param_kappa^2;

vec_x = nominal_x;
diff_x = Inf;
Iter = 0;

fprintf(' ==  Robust RBDO optimization >>> \n');

while (diff_x > param.terminate) && (Iter < 100)
    matK  = sqrtK * diag(vec_x) * sqrtK';
    vec_u = matK \ vec_p;
    grad_comp = -(sqrtK' * vec_u) .* (sqrtK' * vec_u);
    prev_x = vec_x;
    
    cvx_begin sdp quiet
        cvx_solver sedumi
        cvx_precision best
        variable var_z(1,1)
        variable vec_x(nm,1)
        variable g_val(1,1)
        variable mat_W(nm,nm) symmetric
        minimize( dll' * vec_x )
        subject to
            vec_x >= 0;
            g_val + (grad_comp' * tilde_mu)...
                + (param.alpha * norm(grad_comp,1))...
                + trace(mat_W * tilde_Sigma)...
                + (param.beta * norm(reshape(mat_W,nm^2,1),1))...
                + (param_kappa_sq * var_z) <= max_compliance;
            [g_val, vec_p';
                vec_p, (sqrtK * diag(vec_x) * sqrtK')] >= 0;
           [mat_W, grad_comp/2;
                grad_comp'/2, var_z] >= 0;
    cvx_end
    
    Iter = Iter + 1;
    diff_x = norm(vec_x - prev_x);
    
    fprintf('  Iter %g:  obj=%6.5f   change=%3.5e\n',...
        Iter, dll' * vec_x, diff_x);
end
fprintf('     --- \n');
matK = sqrtK * diag(vec_x) * sqrtK';
% <---
% Robust RBDO  optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the results
% --->
fprintf('     uncertainty with L_inf-norm model \n');
fprintf('         alpha = %3.6f ; beta = %3.6f \n',...
    param.alpha, param.beta );
fprintf('         reliability_epsilon = %3.6f \n',...
    param.eps );
fprintf('         termination_threshold = %3.5e \n',...
    param.terminate );
fprintf('     compliance = %6.5f [J] ; bound = %6.2f [J} ; ||p||_inf = %3.3f [kN] \n',...
    (vec_p' * (matK \ vec_p) * 10) , (max_compliance * 10), norm(vec_p,inf) * 10);
fprintf('     opt.val. = %3.5e [cm^3] = %3.5e [mm^3] \n',...
    (dll' * vec_x * 100), (dll' * vec_x * (10^5)) );
fprintf('     opt.sol. = (%6.5f, %6.5f) [cm^2] \n',...
    vec_x(1), vec_x(2) );
fprintf('     nom.sol. = (%6.5f, %6.5f) [cm^2] \n',...
    nominal_x(1), nominal_x(2) );
% <---
% Output the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verification: statistics
% --->
matK  = sqrtK * diag(vec_x) * sqrtK';
vec_u = matK \ vec_p;
g_val = vec_p' * vec_u;
grad_comp = -(sqrtK' * vec_u) .* (sqrtK' * vec_u);

rng(30,'twister');
rand_seed1 = randi(10^6,param.num_data,1);
rng(40,'twister');
rand_seed2 = randi(10^6,param.num_data,1);
rng(50,'twister');
rand_seed3 = randi(10^6,param.num_data,1);

sample_opt_negative = zeros(1,param.num_data);
sample_vec_mu    = zeros(nm,param.num_data);
sample_eig_Sigma = zeros(nm,param.num_data);
failure_prob = zeros(1,param.num_data);
failure_prob_exact = zeros(1,param.num_data);

for iSample = 1:param.num_data
    rng(rand_seed1(iSample),'twister');
    vec_z1 = rand(nm,1) - 0.5;
    vec_z1 = param.alpha * (vec_z1 / norm(vec_z1,Inf));
    rng(rand_seed2(iSample),'twister');
    rand_mat = rand(nm,nm);
    rand_mat = (0.5*(rand_mat+rand_mat')) - 0.5;
    rand_mat = param.beta...
        * ( rand_mat / (norm(reshape(rand_mat,nm^2,1),Inf)) );
    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        variable mat_Z2(nm,nm)
        minimize( norm( mat_Z2 - rand_mat, 'fro' ) )
        subject to
            mat_Z2 >= -param.beta;
            mat_Z2 <=  param.beta;
                tilde_Sigma + mat_Z2 == semidefinite(nm);
    cvx_end
    mat_Z2 = full(mat_Z2);
    
    vec_mu    = tilde_mu    + vec_z1;
    mat_Sigma = tilde_Sigma + mat_Z2;
    if min(eig(mat_Sigma)) < 0
        mat_Sigma = tilde_Sigma;
    end
    
    sample_vec_mu(:,iSample)    = vec_mu;
    sample_eig_Sigma(:,iSample) = sort(eig(mat_Sigma));
    
    sample_opt_negative(iSample) =...
        g_val + (grad_comp' * vec_mu)...
        + (param_kappa * sqrt(grad_comp' * mat_Sigma * grad_comp) )...
        - max_compliance;
    
    rng(rand_seed3(iSample),'twister');
    data_vec_zeta = repmat(vec_mu, 1,param.num_data2)...
        + (randn(param.num_data2,nm) * chol(mat_Sigma))';
    
    sample_compliance   = zeros(param.num_data2,1);
    sample_app_compl    = zeros(param.num_data2,1);
    for jSample = 1:param.num_data2
        sample_x = vec_x + data_vec_zeta(:,jSample);
        sample_K  = sqrtK * diag(sample_x) * sqrtK';
        sample_u = sample_K \ vec_p;
        sample_g_val = vec_p' * sample_u;
        sample_compliance(jSample) = sample_g_val;
        sample_app_compl(jSample) = ...
            g_val + (grad_comp' * data_vec_zeta(:,jSample));
    end
    failure_prob(iSample) =...
        length(find(sample_app_compl > max_compliance)) / param.num_data2;
    failure_prob_exact(iSample) =...
        length(find(sample_compliance > max_compliance)) / param.num_data2;
end
fprintf('     max. L.H.S. of reliability constraint = %3.6f <= 0 \n',...
    max(sample_opt_negative) );
fprintf('     max. of prob. of failure = %3.6f \n',...
    max(failure_prob) );
fprintf('         w/o approximation    = %3.6f \n',...
    max(failure_prob_exact) );
fprintf('     #moment sample = %g ; #zeta sample = %g \n',...
    param.num_data, param.num_data2 );
fprintf(' =========================================================== \n');
% <---
% Verification: statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw figures
% --->
figure;
histogram( failure_prob );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('$\mathrm{P}\{g(\tilde{\mbox{\boldmath $x$}}) + \nabla g(\tilde{\mbox{\boldmath $x$}})^{\top}\mbox{\boldmath $\zeta$}>0\}$', 'Interpreter', 'latex');
xlim([-0.0005, 0.0105]);

figure;
histogram( sample_app_compl - max_compliance);
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('$g(\tilde{\mbox{\boldmath $x$}}) + \nabla g(\tilde{\mbox{\boldmath $x$}})^{\top}\mbox{\boldmath $\zeta$}$  for a distribution', 'Interpreter', 'latex');

figure;
histogram( failure_prob_exact );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('$\mathrm{P}\{g(\tilde{\mbox{\boldmath $x$}}+\mbox{\boldmath $\zeta$})>0\}$', 'Interpreter', 'latex');

figure;
histogram( sample_vec_mu );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Mean', 'Interpreter', 'latex');

figure;
histogram( sample_eig_Sigma );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Eigenvalues of var.-covar.\ matrix', 'Interpreter', 'latex');

figure;
histogram( vec_x(1) + data_vec_zeta(1,:) );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('$x_{1}$ for a distribution', 'Interpreter', 'latex');

figure;
histogram( vec_x(2) + data_vec_zeta(2,:) );
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('$x_{2}$ for a distribution', 'Interpreter', 'latex');
% <---
% Draw figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
