function [lhs, relax_term, imbalance] = energy_balance_diagnostics(psi, n, theta, Jn, Sn, q, VT, tau_E, h)
    % LHS integral
    E = -(psi(2:end) - psi(1:end-1)) / h;
    JnE = Jn .* E;

    lhs = (JnE(1)/2+ sum(JnE(2:end-1)) + JnE(end)/2)*h;

    %RHS 
    relax_term = 1.5*q*n(2:end-1).*(theta(2:end-1) - VT)/tau_E;
    
    int_relax_term = h*(relax_term(1)/2 + ...
        sum(relax_term(2:end-1)) + relax_term(end)/2); 
    
    rhs = Sn(end) - Sn(1) + int_relax_term;

    %Identity check: eps is added in case LHS = 0 as for Va = 0 
    imbalance = abs(lhs - rhs) / (abs(lhs) + eps);