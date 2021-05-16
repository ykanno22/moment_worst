%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2021 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
max_compliance = 100.0;
param.load = 10.0;
eeee = 2.0;
%
param.terminate = 10^(-6);
param.min_cs =  2.0;
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data of the structure
% --->
[dll,matH,coord_x,ir,irr,ird] = member(3,2,1.0);
%
nk = size(coord_x,1); num.node   = nk;
nd = size(matH,1);    num.degree = nd;
nm = size(matH,2);    num.member = nm;
% <---
% Data of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
Idx.deg_of_load = [4, 6];
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
        nominal_x >= param.min_cs;
        [max_compliance, vec_p';
            vec_p, (sqrtK * diag(nominal_x) * sqrtK')] >= 0;
cvx_end

nom_matK = sqrtK * diag(nominal_x) * sqrtK';
fprintf('     compliance = %6.5f [J] ; bound = %6.5f [J} \n',...
    (vec_p' * (nom_matK \ vec_p) * 10) , (max_compliance * 10) );
fprintf('     opt.val. = %3.5e [cm^3] = %3.5e [mm^3] \n',...
    (dll' * nominal_x * 100), (dll' * nominal_x * (10^5)) );

[~] = draw_cs_specified_width(coord_x, irr, nominal_x);
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
        minimize( (dll' * vec_x)...
            +  (1.0 * norm(vec_x - prev_x,2)) )
        subject to
            vec_x >= param.min_cs;
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
fprintf('     nom.val. = %3.5e [cm^3] = %3.5e [mm^3] \n',...
    (dll' * nominal_x * 100), (dll' * nominal_x * (10^5)) );
fprintf('     nom.sol. = (%6.5f, %6.5f) [cm^2] \n',...
    nominal_x(1), nominal_x(2) );
fprintf(' =========================================================== \n');
% <---
% Output the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw figures
% --->
[~] = draw_cs_specified_width(coord_x, irr, vec_x);
% <---
% Draw figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
