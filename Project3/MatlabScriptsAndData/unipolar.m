clear; clc;% close all;

% Physical parameters and mesh
data;

%% Poisson matrix for -eps * psi''
%Nint is used because boundary have already been set
e = ones(Nint,1);
Apoisson = epsSi*spdiags([-e,2*e,-e],-1:1,Nint,Nint)/h^2;

%% Solver parameters
tol       = 1e-10;
maxit     = 500;

%relaxation parameters
omega_n   = 1;
omega_psi = 1;
omega_theta = 1; 

%% Storage

solutions = struct();
solutions_iso_Assignment = struct();
%% Step 1: Initial guess at first bias
Va = Va_list(1);

% Load of the isotermal baseline which are going to be used as initial guess
if isfile('isothermal_baseline.mat')
    load('isothermal_baseline.mat', 'solutions_iso');
    disp('Initial guess from isothermal solution correcly loaded!');
else
    error('isothermal_baseline.mat not found. Run isothermal case first!');
end

if isfile('isothermal_baseline_Assignment.mat')
    load('isothermal_baseline_Assignment.mat', 'solutions_iso_Assignment');
    disp('Baseline from isothermal solution correcly loaded!');
else
    error('isothermal_baseline.mat not found. Run isothermal case first!');
end



% Boundary Conditions for the potential
psi_left  = 0;
psi_right = Va;

% Initial guess for the electron thermal voltage
theta = VT * ones(size(x));

% Boundary condition for the electron carrier density
n(1)   = n_left;
n(end) = n_right;

%% Bias continuation loop:
% Here the initialization can be choosen:
% Set init_strategy = 'I' (isothermal) or 'C' (continuation)
init_strategy = 'C';

for ibias = 1:length(Va_list)
    Va = Va_list(ibias);
    psi_left  = 0;
    psi_right = Va;
    fprintf('\n=========================================\n');
    fprintf('Solving for Va = %.6g V\n', Va);
    fprintf('=========================================\n');
    if init_strategy == 'I'
        %% Case I: Use of isothemal solutions for psi and n
        % Then bias continuation for theta    
        % Initialization with the isotermal solution
        psi = solutions_iso(ibias).psi(:);
        psi(1)   = psi_left;
        psi(end) = psi_right;
    
        n = solutions_iso(ibias).n(:);
    
        Jn = solutions_iso(ibias).Jn(:); %this step is useless: just to be safe...
    
        if ibias == 1
            %initialized before
        else 
            % First guess updated with the previous applied bias
            theta = solutions(ibias-1).theta;
        end
    
    elseif init_strategy == 'C'
        %% Case C: use of the previous bias loop as initial guess of the following one
        if ibias == 1
            psi = VT*log(ND/ND_high) + Va*x/L; 
            
            psi(1)   = psi_left;
            psi(end) = psi_right;
    
            n = ND;
            n(1)   = n_left;
            n(end) = n_right;
            theta = VT * ones(size(x));
        else
            % Use previous solution as initial guess.
            n = solutions(ibias-1).n; 
            n(1)   = n_left;
            n(end) = n_right;

            % Affine correction of potential to match the new voltage.
            Va_old = Va_list(ibias-1);
            dVa = Va - Va_old;
    
            psi = solutions(ibias-1).psi + dVa*x/L;
            psi(1)   = psi_left;
            psi(end) = psi_right;
    
            theta = solutions(ibias-1).theta;
            theta(1)   = theta_left;
            theta(end) = theta_right;
        end
    end
    err = 1 + tol;
    it  = 0;

    while err > tol && it < maxit

        it = it + 1;

        psi_old = psi;
        n_old   = n;

        %% Step 2: Continuity equation with Scharfetter-Gummel
        % find n
        n_new = continuity_n(psi, theta, mun, n_left, n_right, Nint,h);

        % Damping on density
        n = (1-omega_n)*n_old + omega_n*n_new;

        % Enforce density boundary values
        n(1)   = n_left;
        n(end) = n_right;
        
        %% Step 3: Jn calculation
        Jn = electron_current_sg(psi, n, h, theta, mun, q);
       
        %% Step 4: Solving of energy equation
        % Averages 
        n_avg = (n(1:end-1)+n(2:end))/2;
        theta_avg = (theta(1:end-1) + theta(2:end))/2;
       
        % lambda vector calculation
        lambda = c_lambda*q*mun.*n_avg.*theta_avg;
        
        theta_old = theta;
        
        [theta_new, Sn, H] = energy_theta(psi, n, Jn, lambda, theta_left, theta_right, Nint, h, tau_E, VT, q);

        theta = (1-omega_theta)*theta_old + omega_theta*theta_new;    
        %enforcing BC
        theta(1) = theta_left;
        theta(end) = theta_right;
        
        %If the following line gets uncommentend, isothermal solutions are
        %computed:
        %theta = VT*ones(size(x));
        
        %% 5. Gummel-linearized Poisson equation
        % Applied linearization:
        %
        % n(psi^{k+1}) approx n^k + (n^k/VT)(psi^{k+1}-psi^k)
        %
        % Therefore:
        %
        % -eps psi'' + q n^k/VT psi
        % =
        % q(ND - n^k) + q n^k/VT psi^k
        %From which the matrix have been extracted
        
        n_int      = n(2:end-1);
        ND_int     = ND(2:end-1);
        
        psi_old_in = psi_old(2:end-1);
        
        theta_int = theta(2:end-1);

        Mreact = spdiags(q*n_int./theta_int,0,Nint,Nint);
        
        %RHS
        b = q*(ND_int - n_int) + q*(n_int./theta_int).*psi_old_in;
       
        % Dirichlet boundary contributions
        b(1)   = b(1)   + epsSi*psi_left/h^2;
        b(end) = b(end) + epsSi*psi_right/h^2;
        
        %computation on the new psi
        psi_new = [psi_left; (Apoisson + Mreact)\b; psi_right];

        % Damping on potential
        psi = (1-omega_psi)*psi_old + omega_psi*psi_new;

        % Enforce potential boundary values
        psi(1)   = psi_left;
        psi(end) = psi_right;

        %% Step 7. Error indicators

        err_psi = norm(psi - psi_old,inf)/max([VT,norm(psi,inf)]);
        err_n   = norm(n - n_old,inf)/max(norm(n,inf),1);     
        err_theta = norm(theta - theta_old, inf)/ norm(theta, inf);
        
        err = max([err_psi, err_n, err_theta]);

        fprintf('it = %4d | err = %.3e | err_psi = %.3e | err_n = %.3e | err_th = %.3e | min psi = %.4e | max psi = %.4e\n', ...
            it, err, err_psi, err_n, err_theta, min(psi), max(psi));

    end

    if it == maxit
        warning('Maximum number of iterations reached at Va = %.4g V',Va);
    end

    %% Save solution
    solutions(ibias).Va  = Va;
    solutions(ibias).psi = psi;
    solutions(ibias).n   = n;
    solutions(ibias).Jn  = Jn;
    solutions(ibias).theta = theta;
    
    
    
    %% Plot current solution
    figure(1); clf;

    subplot(6,1,1);
    semilogy(x*1e6,ND,'k--','LineWidth',1.2); hold on;
    semilogy(x*1e6,n,'b','LineWidth',1.4);
    semilogy(x*1e6, solutions_iso_Assignment(ibias).n, 'r','LineWidth',1.4);
    grid on;
    xlabel('x [\mum]');
    ylabel('density [m^{-3}]');
    legend('N_D', 'n (Non-Iso)', 'n (Iso)', 'Location', 'best');
    title(sprintf('Electron density, Va = %.4g V',Va));

    subplot(6,1,2);
    plot(x*1e6,psi,'LineWidth',1.4); hold on;
    plot(x*1e6, solutions_iso_Assignment(ibias).psi, 'r','LineWidth',1.4);
    grid on;
    xlabel('x [\mum]');
    ylabel('\psi [V]');
    legend('\psi (Non-Iso)', '\psi (Iso)', 'Location', 'best');
    title('Electrostatic potential');

    subplot(6,1,3);
    Jn = solutions(ibias).Jn; hold on;
    plot(x(1:end-1)*1e6,Jn,'LineWidth',1.4);
    plot(x(1:end-1)*1e6, solutions_iso_Assignment(ibias).Jn, 'r','LineWidth',1.4);
    grid on;
    xlabel('x [\mum]');
    ylabel('J_n [A/m^2]');
    legend('J_n (Non-Iso)', 'J_n (Iso)', 'Location', 'best');
    title('SG electron current density');

    %Plot for electron temperature
    Tn = theta * q / kB;

    subplot(6,1,4)
    plot(x*1e6, Tn, 'b', 'LineWidth', 1.4);
    hold on;
    yline(T, 'k--', 'T_L (Ambient)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    grid on;
    xlabel('x [\mum]');
    ylabel('T_n [K]');
    title(sprintf('Electron Temperature, Va = %.3g V', Va));

    % % plot for theta    
    % subplot(6,1,4)
    % plot(x*1e6,theta,'r','LineWidth',1.4);
    % yline(VT, 'k--', 'V_T (Ambiente)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    % grid on;
    % hold on;
    % xlabel('x [\mum]');
    % ylabel('\vartheta [V]');
    % title(sprintf('Electron Thermal Voltage, Va = %.3g V',Va));

    subplot(6,1,5)
    plot(x*1e6,n./ND,'LineWidth',1.4); hold on;
    plot(x*1e6,(solutions_iso_Assignment(ibias).n)./ND,'LineWidth',1.4);
    grid on;
    xlabel('x [\mum]');
    ylabel('n/N_D');
    legend('n/N_D (ET)', 'n/N_D (Iso)', 'Location', 'best');
    title(sprintf('Quasi-neutrality ratio, Va = %.3g V',Va));

    subplot(6,1,6)
    plot(x*1e6,q*(ND - n),'LineWidth',1.4); hold on;
    plot(x*1e6,q*(ND - (solutions_iso_Assignment(ibias).n)),'LineWidth',1.4);   
    grid on;
    xlabel('x [\mum]');
    ylabel('\rho = q(N_D-n) [C/m^3]');
    legend('\rho (Non-Iso)', '\rho (Iso)', 'Location', 'best');
    title(sprintf('Space charge density, Va = %.3g V',Va));

    drawnow;
    %% Final residual diagnostics
    fprintf('\nFinal diagnostics for Va = %.6e V\n', Va);

    fprintf('min(psi)                  = %.6e V\n', min(psi));
    fprintf('max(psi)                  = %.6e V\n', max(psi));

    fprintf('min(n)                    = %.6e m^-3\n', min(n));
    fprintf('max(n)                    = %.6e m^-3\n', max(n));
    fprintf('min(n/ND)                 = %.6e\n', min(n./ND));
    fprintf('max(n/ND)                 = %.6e\n', max(n./ND));

    fprintf('mean(Jn)                  = %.6e A/m^2\n', mean(Jn));
    fprintf('min(Jn)                   = %.6e A/m^2\n', min(Jn));
    fprintf('max(Jn)                   = %.6e A/m^2\n', max(Jn));
    
    % --------- Current conservation diagnostic -----------%
    %Jn is computed again because at the end of the Gummel iterations the
    %last theta have to be taken into account
    Jn = electron_current_sg(psi, n, h, theta, mun, q);
    Jspread = max(Jn) - min(Jn);
    
    Jscale = q*mun*VT*ND_high/L;
    
    rel_scale   = Jspread / Jscale;
    
    fprintf('\n--- Diagnostic 11.1: Current conservation ---\n');
    fprintf('Jscale                        = %.6e A/m^2\n', Jscale);
    fprintf('spread(Jn)                    = %.6e A/m^2\n', Jspread);
    fprintf('spread(Jn) / Jscale           = %.6e\n', rel_scale);
    
    % This ratio is meaning full only far form the equilibrium
    if abs(mean(Jn)) > 1e-12 * Jscale
        rel_mean = Jspread / abs(mean(Jn));
        fprintf('spread(Jn) / |mean(Jn)|       = %.6e\n', rel_mean);
    else
        fprintf('spread(Jn) / |mean(Jn)|       = not significative, current is almost nil\n');
    end
    
    %-------------- Poisson residual Diagnostic ---------------%

    rho_int = q*(ND(2:end-1) - n(2:end-1));

    lhsP = -epsSi*(psi(3:end)-2*psi(2:end-1)+psi(1:end-2))/h^2;
    rhsP = rho_int;

    resP = lhsP - rhsP;

    Pscale = q*ND_high;
    fprintf('\n--- Diagnostics 11.2: Poisson Residual ---\n');
    fprintf('Poisson residual inf      = %.6e C/m^3\n', norm(resP,inf));
    fprintf('Poisson residual / Pscale = %.6e\n', norm(resP,inf)/Pscale);
    fprintf('Poisson residual / max(abs(rho_int)) = %.6e\n', ...
        norm(resP,inf)/max(norm(rho_int,inf),1));
    
    %----------------- Energy residual Diagnostic -------------%
    
    alpha = 1.5*q*n(2:end-1)/tau_E;
    
    dSn = (Sn(2:end)-Sn(1:end-1))/h;
    
    resEN = dSn - (H - alpha.*(theta_int - VT));
    
    %denom_EN_RES = max(abs(H) + abs(alpha.*(theta_int - VT))) + eps; 
    denom_EN_RES = max(abs(H + alpha.*(theta_int - VT))) + eps; 

    fprintf('\n--- Diagnostics 11.3: Energy residual ---\n')
    fprintf('||R^E||_inf / denom           = %.6e\n', norm(resEN, inf) / denom_EN_RES);
    
    % ---------------- Global Energy Balance Diagnostic ----------------%

    [~, ~, imbalance] = energy_balance_diagnostics(psi, n, theta, Jn, Sn, q, VT, tau_E, h);
    
    fprintf('\n--- Diagnostics 11.4: Global energy balance ---\n');
    fprintf('Relative imbalance            = %.6e\n', imbalance);
    
    %----------------- Drift Velocity Diagnostics --------------%
   
    vsat = 1e5;   % m/s, typical electron saturation velocity in Si
    n_edge = 0.5*(n(1:end-1) + n(2:end));

    % Actual carrier mean velocity
    u_edge = -Jn./(q*n_edge);
    fprintf('\n--- Diagnostics 11.5: Velocity and Validity indicators ---\n');
    
    fprintf('max(|u_n|)              = %.6e m/s\n', max(abs(u_edge)));
    fprintf('max(|u_n|)/vsat         = %.6e\n', max(abs(u_edge))/vsat);
    
    ratio_vsat = max(abs(u_edge))/vsat;
    if ratio_vsat < 0.2
        fprintf('Regime: low-field, mobility model valid\n');
    elseif ratio_vsat < 0.3
        fprintf('Regime: acceptable \n');
    else
        fprintf('Regime: WARNING - constant mobility questionable\n');
    end
    
    %----------------- Quasi-Fermi diagnostics --------------%
    %it is applied for the isothermal case: may be useful to compare with 
    % velocity drift checks
    psi_iso = solutions_iso(ibias).psi;
    n_iso   = solutions_iso(ibias).n;

    psiF  = psi_iso - VT * log(n_iso / ND_high);
    Fqf   = -diff(psiF) / h;
    u_qf  = mun * abs(Fqf);

    fprintf('max(mu*|grad psiF|)/vsat = %.6e  [isothermic]\n', max(u_qf)/vsat);    

    %---------------Temperature Checks -------------------------%
    
    % Electron temperature
    Tn = theta * q / kB;

    Tn_max = max(Tn);
    Tn_min = min(Tn);
    dT_rel_max = max((Tn-T)/T);
    
    fprintf('\n--- Diagnostics 11.6: Temperature checks ---\n');
    fprintf('max(Tn)                       = %.4f K\n', Tn_max);
    fprintf('min(Tn)                       = %.4f K\n', Tn_min);
    fprintf('max((Tn-TL)/TL)               = %.6e\n',  dT_rel_max);

    if Tn_min <= 0
        error('The temperature is negative!');
    end
    %keyboard;
end
%-------------- Report on carrier velocity for every bias-----------------%
fprintf('\n--- Report on carrier velocity for every bias ---\n');
fprintf('%-10s %-15s %-15s\n', 'Va [V]', 'max|u|/vsat', 'regime');
for ibias = 1:length(solutions)
    Va   = solutions(ibias).Va;
    Jn   = solutions(ibias).Jn;
    n    = solutions(ibias).n;
    
    n_edge = 0.5*(n(1:end-1)+n(2:end));
    u_edge = -Jn./(q*n_edge);
    
    ratio  = max(abs(u_edge))/vsat;
    if ratio < 0.2,      regime = 'low-field';
    elseif ratio < 0.3,  regime = 'acceptable';
    else,                regime = 'WARNING';
    end
    fprintf('%-10.3f %-15.3e %-15s\n', Va, ratio, regime);
end

%% I-V characteristics
Va_values = arrayfun(@(s) s.Va, solutions);
J_values  = arrayfun(@(s) mean(s.Jn), solutions);

J_iso_values = arrayfun(@(s) mean(s.Jn), solutions_iso_Assignment);

figure(2); clf;
plot(Va_values,J_values, 'o-r','LineWidth',1.4);
hold on;
plot(Va_values,J_iso_values,'o-b','LineWidth',1.4);
grid on;
xlabel('V_a [V]');
ylabel('mean J_n [A/m^2]');
legend('non_isothermal', 'Isothermal', 'Location', 'best');
title('Current-voltage characteristic');

%% Current comparison: Non-Isothermal vs Isotermal (Va = 0 excluded)
% Length extraction 
N_sim = length(solutions);


if N_sim > 1
    %I ignore Va=0
    Va_values     = arrayfun(@(s) s.Va, solutions(2:N_sim));
    J_values      = arrayfun(@(s) mean(s.Jn), solutions(2:N_sim));
    J_iso_values = arrayfun(@(s) mean(s.Jn), solutions_iso_Assignment(2:N_sim));

    % absolute difference computation
    Delta_J = J_values - J_iso_values; 

    % relative difference computation
    Delta_J_perc = zeros(1, length(Va_values));
    idx_nz = abs(J_iso_values) > 1e-12; 
    Delta_J_perc(idx_nz) = (Delta_J(idx_nz) ./ abs(J_iso_values(idx_nz))) * 100;

    %Plot
    figure(3); clf;
    set(gcf, 'Name', 'Comparison ET vs DD');

    % Subplot 1: absolute difference
    subplot(2,1,1);
    plot(Va_values, Delta_J, 'o-', 'LineWidth', 1.5, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
    grid on;
    xlabel('V_a [V]');
    ylabel('\Delta J_n [A/m^2]');
    title('Absolute Difference: J_{n, ET} - J_{n, DD} (V_a > 0)');

    % Subplot 2: relative difference
    subplot(2,1,2);
    plot(Va_values, Delta_J_perc, 's-', 'LineWidth', 1.5, 'Color', '#0072BD', 'MarkerFaceColor', '#0072BD');
    grid on;
    xlabel('V_a [V]');
    ylabel('\Delta J_n / J_{n, DD} [%]');
    title('Relative variation of current (V_a > 0)');
 else
     disp('Warning: only Va = 0 have been computed. No comparison available');
 end

%% Figure 4: Mesh refinement (persistent across runs, single bias)
ibias = 7;
Va_target = Va_list(ibias);
Jn = solutions(ibias).Jn;

figure(4);
hold on;
plot(x(1:end-1)*1e6, Jn, 'LineWidth', 1.4, ...
     'DisplayName', sprintf('N = %d', N));
grid on;
xlabel('x [\mum]');
ylabel('J_n [A/m^2]');
title(sprintf('Mesh refinement: J_n(x),  V_a = %.2g V', Va_target));
legend('Location', 'best');
drawnow;

ibias = 4;
Va_target = Va_list(ibias);
Jn = solutions(ibias).Jn;

figure(5);
hold on;
plot(x(1:end-1)*1e6, Jn, 'LineWidth', 1.4, ...
     'DisplayName', sprintf('N = %d', N));
grid on;
xlabel('x [\mum]');
ylabel('J_n [A/m^2]');
title(sprintf('Mesh refinement: J_n(x),  V_a = %.2g V', Va_target));
legend('Location', 'best');
drawnow;
