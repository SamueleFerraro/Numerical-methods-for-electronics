function B = BernoulliF(z)

B = zeros(size(z));

small = abs(z) < 1e-8;
large = ~small;

% Series expansion near zero:
% B(z) = z/(exp(z)-1)
%      = 1 - z/2 + z^2/12 - z^4/720 + ...
zs = z(small);
B(small) = 1 - zs/2 + zs.^2/12 - zs.^4/720;

zl = z(large);
B(large) = zl ./ expm1(zl);

end