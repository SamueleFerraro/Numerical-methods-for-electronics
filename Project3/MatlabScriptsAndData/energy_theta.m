function [theta, Sn, H] = energy_theta(psi, n, Jn, lambda, theta_a, theta_b, Nint, h, tau_E, VT, q)

% Peclet number

PE = (5*Jn*h)./(2*lambda);
[bp,bn] = bern(PE);

% %electric field computation
 E = -(psi(2:end)-psi(1:end-1))/h;
% E1 = -diff(psi)/h;


% H coefficient of the equation 
H = 0.5*(Jn(2:end).*E(2:end) + Jn(1:end-1).*E(1:end-1));

% alpha coefficient: only internal node of n are used for alpha vector
alpha = 1.5*(q*n(2:end-1))/tau_E;

diag_main_en = (lambda(2:end).*bn(2:end)+lambda(1:end-1).*bp(1:end-1))/h;
diag_main = diag_main_en + h*alpha;

% Subdiagonal: row r couples to node r-1
% Valid for r = 2,...,Nint
rows_low = (2:Nint).';
cols_low = (1:Nint-1).';
vals_low = -(lambda(2:Nint).*bn(2:Nint))/h ;

% Superdiagonal: row r couples to node r+1
% Valid for r = 1,...,Nint-1
rows_up = (1:Nint-1).';
cols_up = (2:Nint).';
vals_up = -(lambda(2:Nint).*bp(2:Nint))/h;

% Main diagonal
rows_diag = (1:Nint).';
cols_diag = (1:Nint).';
vals_diag = diag_main(:);

% Assemble matrix
rows = [rows_diag; rows_low; rows_up];
cols = [cols_diag; cols_low; cols_up];
vals = [vals_diag; vals_low; vals_up];

A = sparse(rows,cols,vals,Nint,Nint);

%RHS of the equation
rhs = h*(H+alpha*VT);


A_low_boundary = -(lambda(1)*bn(1))/h;
A_up_boundary  = -(lambda(end)*bp(end))/h;

% I must include external nodes (Boundary conditions)
rhs(1)   = rhs(1)   - A_low_boundary * theta_a;
rhs(end) = rhs(end) - A_up_boundary  * theta_b;

% Solve the linear system for theta
theta_int = A \ rhs; 

theta = [theta_a; theta_int; theta_b];

%%For diagnostic 11.3
Sn = (lambda)/h .* (theta(1:end-1).*bn-theta(2:end).*bp); 

if any(theta < 0)
        error('!!! Negative thetha !!!');
end


