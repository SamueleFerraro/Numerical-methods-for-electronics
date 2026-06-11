function n = continuity_n(psi, theta, mun, na, nb, Nint, h)
%% Used to solve current continuity equation

% Nint = number of interior nodes
% total number of nodes is Nint+2

psi   = psi(:);
theta = theta(:);
%Computation of the logarithmic mean of theta: 
theta_log_mean = log_mean(theta(2:end), theta(1:end-1));

% theta = diff(psi)/Vth;

eta = (diff(psi) - diff(theta)) ./ theta_log_mean;

% Convention:
% bp(i) = B(eta_i)
% bn(i) = B(-eta_i)
[bp,bn] = bern(eta);

Dn = mun*theta_log_mean;
%c is a vector in the non-isothermal case
c = Dn/h;

% Diagonal entries for interior nodes
diag_main = c(1:end-1).*bp(1:end-1) + c(2:end).*bn(2:end);

% Subdiagonal: row r couples to node r-1
% Valid for r = 2,...,Nint
rows_low = (2:Nint).';
cols_low = (1:Nint-1).';
vals_low = -c(2:end-1).*bn(2:end-1);

% Superdiagonal: row r couples to node r+1
% Valid for r = 1,...,Nint-1
rows_up = (1:Nint-1).';
cols_up = (2:Nint).';
vals_up = -c(2:end-1).*bp(2:end-1);

% Main diagonal
rows_diag = (1:Nint).';
cols_diag = (1:Nint).';
vals_diag = diag_main(:);

% Assemble matrix
rows = [rows_diag; rows_low; rows_up];
cols = [cols_diag; cols_low; cols_up];
vals = [vals_diag; vals_low; vals_up];

An = sparse(rows,cols,vals,Nint,Nint);

% Boundary contributions
Fn = sparse(Nint,1);
Fn(1)   = c(1)*bn(1)*na;
Fn(end) = c(end)*bp(end)*nb;

n_int = An\Fn;

n = [na; n_int; nb];

end