function [p, C] = convergence_rate(err, h)
% this function estimates the order of convergence p and the asymptotic 
% error constant C of a sequence of errors contained in the vector err
% [p, C] = order_estimate(err)
%
p = log2(err(1:end-1)./err(2:end))./log2(h(1:end-1)./h(2:end));
%
C = err(1:end-1)./(h(1:end-1).^p);
return