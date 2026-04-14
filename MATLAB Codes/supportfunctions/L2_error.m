function L2_error = L2_error(x, uh, uex, r)
% Authors: A. Tonini
%
% err = L2_error(x, uh, uex, r)
%
% This script evaluates the L2 error between a FEM approximate 
% solution uh (vector of dofs) evaluated at nodes x and the exact 
% solution uex (function handle)
%
% x: nodes vector
% uh: nodal FEM solution
% uex: exact solution function handle
% r: polynomial degree

% Number of elements
Ne = (length(x) - 1)/max(r,1);
L2_error = 0;

% Gauss quadratures on a reference element [-1,1].
% You need to transform the gauss points from the reference
% element to the one on which we are integrating
switch r
    case 0
        gauss_pts  = 0;
        gauss_w    = 2;
        % Basis functions on the reference element [-1,1]
        Nphi = @(x) 1.0+0.0*x;
    case 1
        gauss_pts  = [-1/sqrt(3), 1/sqrt(3)];
        gauss_w    = [1, 1]';
        % Basis functions on the reference element [-1,1]
        Nphi = @(x) [0.5*(1-x); 0.5*(1+x)];
    case 2
        gauss_pts  = [-sqrt(3/5), 0, sqrt(3/5)];
        gauss_w    = [5/9, 8/9, 5/9]';
        % Basis functions on the reference element [-1,1]
        Nphi = @(x) [0.5*x.*(x-1); 1-x.^2; 0.5*x.*(x+1)];
    otherwise
        error("FEM for polynomial degree %d not implemented.\n", r);
end
Nphi_q = Nphi(gauss_pts);

Nq = length(gauss_w);
local_dofs = r+1;
for e = 1:Ne
    if r == 0
        x1 = x(e);
        x2 = x(1+e);

        % Local degrees of freedom
        uh_local = uh(e);
    else
        x1 = x(1+r*(e-1));
        x2 = x(1+r*e);
        
        % Local degrees of freedom
        uh_local = uh(1+r*(e-1):1+r*e);    
    end
    he = x2 - x1;
    
    % Mapping from the reference element [-1,1] to the element
    % [x(1+r*(e-1)), x(1+r*e)]
    xq = 0.5*(x2-x1)*gauss_pts+0.5*(x2+x1);
    
    % Jacobian of the transformation for the integral on the element e
    J  = he/2;
        
    % Cycle over the quadrature nodes
    for q = 1:Nq
        u_q_local = 0;
         % Cycle over the basis functions
        for i = 1:local_dofs
            u_q_local = u_q_local + uh_local(i)*Nphi_q(i,q);
        end
        L2_error = L2_error + J * gauss_w(q) * (u_q_local - uex(xq(q)))^2;
    end
end

L2_error = sqrt(L2_error);

end