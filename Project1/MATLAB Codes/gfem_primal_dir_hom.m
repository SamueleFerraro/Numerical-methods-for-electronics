function [x, u, J] = gfem_primal_dir_hom(a, b, N, f, mu, beta, sigma)
% [x, u, J] = gfem_primal_dir(a, b, N, f, mu, beta, sigma)
%
% This script implements GFEM for the linear advection-diffusion-reaction
% boundary problem:
%
% -(mu(x) u'(x))' + mu u'(x) + sigma*u(x) = f(x)    x in (a,b)
% u(a) = 0
% u(b) = 0
% 
% J is computed in post processing as formulated in Point a) of the report

% space discretization
h = (b-a)/(N+1);
x = linspace(a, b, N+2)';
u = zeros(N+2, 1);

% Mean points of x
xm = 0.5*(x(1:end-1)+x(2:end));

% assembly
e = ones(N,1);

%computation of mu mean via trapezoidal quadrature formula
mu_mean = 0.5*(mu(x(1:end-1)) + mu(x(2:end)));

%Main diagonal of the stiffness matrix: calculated at lecture 2 notes
mu_diag = mu_mean(1:end-1) + mu_mean(2:end);

K = 1/h * spdiags([-[mu_mean(2:end-1); 0], mu_diag, -[0; mu_mean(2:end-1)]], [-1, 0, 1], N, N);

B = beta*0.5 * spdiags([-e, e], [-1, 1], N, N);
M = sigma*h/6 * spdiags([e, 4*e, e], [ -1, 0, 1], N, N);

% We can avoid the modification of the system matrices on the boundary nodes 
% as we are applying Dirichlet BCs
A = K + B + M;

% Midpoint quadrature rule again with
%Computation the rhs
F = h/2 * f(xm(1:end-1)) + h/2 * f(xm(2:end));


% Solve linear system
u(2:end-1) = A\F;

%Computation post processing of J primal
J = -mu_mean .* (u(2:end)-u(1:end-1))/h;

return
