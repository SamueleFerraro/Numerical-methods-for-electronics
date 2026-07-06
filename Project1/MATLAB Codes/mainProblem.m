clear variables
close all
clc

addpath("./supportfunctions");

%parameters
a = 0;
b = 1;
sigma = 5;
N = 29;

beta=0;

%Functions given by the assignment and the ones needed to build f
mu = @(x) 1.5 + sin(30*pi*x);
dmu = @(x) 30*pi*cos(30*pi*x);

u_ex = @(x) x.^(5/2) .* (1-x);
du_ex = @(x) 2.5*x.^(1.5).*(1-x) - x.^(5/2);
d2u_ex = @(x) 3.75*x.^(0.5).*(1-x) - 5*x.^(1.5);

J_ex    = @(x) -mu(x) .* du_ex(x);

f = @(x) -dmu(x).*du_ex(x) - mu(x).*d2u_ex(x) + sigma*u_ex(x);

%solving the problem using the two different formulations:
[x_prim, u_prim, J_prim] = gfem_primal_dir_hom(a, b, N, f, mu, beta, sigma);
[x_mix, u_mix, J_mix] = gfem_mixed_dir_hom(a, b, N, f,mu , sigma);

xex = linspace(a, b, 1000)';

% --- Plot u(x) ---
setfonts;
figure;
plot(x_mix, u_mix, '-xb', x_prim, u_prim, '-r', xex, u_ex(xex), '--k');
title('Density u(x)');
xlabel('x');
ylabel('u');
legend('Mixed', 'Primal', 'Exact');


% --- Plot J(x) ---
setfonts;
figure;
hold on
for i = 1:length(J_mix)
    plot([x_mix(i) x_mix(i+1)], [J_mix(i) J_mix(i)], "-b", "HandleVisibility","off");
    plot([x_prim(i) x_prim(i+1)], [J_prim(i) J_prim(i)], "-r", "HandleVisibility","off");   
end
xm = 0.5*(x_mix(1:end-1) + x_mix(2:end));
plot(xm, J_mix, 'xb', xm, J_prim, 'or', xex, J_ex(xex), '--k');
title('Flux J(x)');
xlabel('x');
ylabel('J');
legend('Mixed','Primal', 'Exact');

%%Errors
%Calculation of the errors H1 of u
[err_H1_mix, err_semi_H1_mix] = H1_error(x_mix, u_mix, u_ex, du_ex, 1);
[err_H1_prim, err_semi_H1_prim] = H1_error(x_prim, u_prim, u_ex, du_ex, 1);

fprintf('H1 error on u  - Mixed:  %.6e  (seminorm: %.6e)\n', err_H1_mix, err_semi_H1_mix);
fprintf('H1 error on u  - Primal: %.6e  (seminorm: %.6e)\n', err_H1_prim, err_semi_H1_prim);

%Calculation of the errors L2 of J
err_L2_mix = L2_error(xm, J_mix, J_ex, 0);
err_L2_prim = L2_error(xm, J_prim, J_ex, 0);

fprintf('L2 error on J  - Mixed:  %.6e\n', err_L2_mix);
fprintf('L2 error on J  - Primal: %.6e\n', err_L2_prim);

%% Point e) - Convergence analysis

N_vec = [3^5, 3^6, 3^7, 3^8, 3^9, 3^10];
h_vec = ((b-a) ./ (N_vec + 1))';

%Initialization of the Mixed formulation
err_combined_mix = zeros(length(N_vec), 1);
err_H1_u_mix     = zeros(length(N_vec), 1);
err_semi_H1_u_mix = zeros(length(N_vec), 1);
err_L2_J_mix     = zeros(length(N_vec), 1);

%Initialization of the Primal formulation
err_H1_u_prim     = zeros(length(N_vec), 1);
err_semi_H1_u_prim = zeros(length(N_vec), 1);
err_L2_J_prim     = zeros(length(N_vec), 1);

for k = 1:length(N_vec)
    N_k = N_vec(k);

    %Solve mixed formulation
    [x_mix_k, u_mix_k, J_mix_k] = gfem_mixed_dir_hom(a, b, N_k, f, mu, sigma);
    
    % Compute xm for
    xm_mix_k = 0.5*(x_mix_k(1:end-1) + x_mix_k(2:end));

    % Compute xm for J
  
    % H1 error on u (Mixed)
    [err_H1_u_mix(k), err_semi_H1_u_mix(k)] = H1_error(x_mix_k, u_mix_k, u_ex, du_ex, 1);

    % L2 error on J (Mixed)
    err_L2_J_mix(k) = L2_error(xm_mix_k, J_mix_k, J_ex, 0);

    % Combined error (left hand side of theorem)
    err_combined_mix(k) = err_H1_u_mix(k) + err_L2_J_mix(k);
    
    fprintf('\n--- Combined error per iteration ---\n');
    fprintf('N+1 = %6d | h = %.2e | Err_comb = %.6e\n', N_k+1, h_vec(k), err_combined_mix(k));

    % --- Solve primal formulation ---
    [x_prim_k, u_prim_k, J_prim_k] = gfem_primal_dir_hom(a, b, N_k, f, mu, beta, sigma);
    xm_prim_k = 0.5*(x_prim_k(1:end-1) + x_prim_k(2:end));

    % H1 error on u (Primal)
    [err_H1_u_prim(k), err_semi_H1_u_prim(k)] = H1_error(x_prim_k, u_prim_k, u_ex, du_ex, 1);
    % L2 error on J (Primal)
    err_L2_J_prim(k) = L2_error(xm_prim_k, J_prim_k, J_ex, 0);
end

% --- Estimate convergence order (Mixed) ---
[p_H1_u_mix, C_H1_u_mix] = convergence_rate(err_H1_u_mix, h_vec);
[p_L2_J_mix, C_L2_J_mix] = convergence_rate(err_L2_J_mix, h_vec);
[p_combined_mix, C_comb_mix] = convergence_rate(err_combined_mix, h_vec);

% --- Estimate convergence order (Primal) ---
[p_H1_u_prim, C_H1_u_prim] = convergence_rate(err_H1_u_prim, h_vec);
[p_L2_J_prim, C_L2_J_prim] = convergence_rate(err_L2_J_prim, h_vec);

fprintf('\n--- Estimate of convergence order (Mixed) ---\n');
fprintf('H1(u) order:      %.4f. Computed constant C: %.4f \n', p_H1_u_mix(end),C_H1_u_mix(end));
fprintf('L2(J) order:      %.4f. Computed constant C: %.4f \n', p_L2_J_mix(end), C_L2_J_prim(end));
fprintf('Combined order:   %.4f. Computed constant C: %.4f \n', p_combined_mix(end), C_comb_mix(end));

fprintf('\n--- Estimate of convergence order (Primal) ---\n');
fprintf('H1(u) order:      %.4f. Computed constant C: %.4f \n', p_H1_u_prim(end), C_H1_u_prim(end));
fprintf('L2(J) order:      %.4f. Computed constant C: %.4f \n', p_L2_J_prim(end), C_L2_J_prim(end));

%% Plot of the results (Mixed)
setfonts;
figure('Name', 'Convergence Analysis - Mixed Formulation', 'Position', [100, 100, 1600, 500]);

% --- Subplot 1: Combined error (Theorem 1) ---
subplot(1,3,1);
loglog(h_vec, err_combined_mix, 'o-', ...
       h_vec, h_vec, '--', ...
       h_vec, h_vec.^2, '--', ...
        'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Error');
legend('$\|u-u_h\|_V + \|J-J_h\|_Q$', '$p = 1$', '$p = 2$',  ...
       'Interpreter', 'latex', 'Location', 'southeast');
title('Combined Error');

% --- Subplot 2: Error on u (H1 Norm e Seminorm) ---
subplot(1,3,2);
loglog(h_vec, err_H1_u_mix, 'o-', ...
       h_vec, err_semi_H1_u_mix, 's-', ...
       h_vec, h_vec, '--', ...
       h_vec, h_vec.^2, '--', ...
        'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Error');
legend('$\|u-u_h\|_{H^1}$', '$|u-u_h|_{H^1}$', '$p = 1$', '$p = 2$', ...
       'Interpreter', 'latex', 'Location', 'southeast');
title('Error on u');

% --- Subplot 3: Error on J (L2 Norm) ---
subplot(1,3,3);
loglog(h_vec, err_L2_J_mix, 'o-', ...
       h_vec, h_vec, '--', ...
       h_vec, h_vec.^2, '--', ...
         'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Error');
legend('$\|J-J_h\|_{L^2}$', '$p = 1$', '$p = 2$', ...
       'Interpreter', 'latex', 'Location', 'southeast');
title('Error on J');

%% Plot dei risultati - Primal (1x2)
figure('Name', 'Convergence Analysis - Primal Formulation', 'Position', [150, 150, 1100, 500]);

% --- Subplot 1: Errore su u (Primal) ---
subplot(1,2,1);
loglog(h_vec, err_H1_u_prim, 'o-', ...
       h_vec, err_semi_H1_u_prim, 's-', ...
       h_vec, h_vec, '--', ...
       h_vec, h_vec.^2, '--', ...
        'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Error');
legend('$\|u-u_h\|_{H^1}$', '$|u-u_h|_{H^1}$', '$p = 1$', '$p = 2$', ...
       'Interpreter', 'latex', 'Location', 'southeast');
title('Error on u (Primal)');

% --- Subplot 2: Errore su J (Primal) ---
subplot(1,2,2);
loglog(h_vec, err_L2_J_prim, 'o-', ...
       h_vec, h_vec, '--', ...
       h_vec, h_vec.^2, '--', ...
         'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Error');
legend('$\|J-J_h\|_{L^2}$', '$p = 1$', '$p = 2$', ...
       'Interpreter', 'latex', 'Location', 'southeast');
title('Error on J (Primal)');


rmpath('./supportfunctions');