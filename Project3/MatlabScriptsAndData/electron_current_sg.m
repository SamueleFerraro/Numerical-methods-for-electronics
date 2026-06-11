function Jn = electron_current_sg(psi, n, h, theta, mun, q)

%Computation of the logarithmic mean of theta: 

theta_log_mean = log_mean(theta(2:end), theta(1:end-1));

eta = (diff(psi) - diff(theta)) ./ theta_log_mean;

Dn = mun*theta_log_mean;

[bp,bn] = bern(eta);

Jn = q*Dn/h.*(n(2:end).*bp - n(1:end-1).*bn );

end