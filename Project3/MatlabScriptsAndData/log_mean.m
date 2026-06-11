function lm = log_mean(a,b)
%log_mean computed the logarithmic mean between two vectors a & b
a=a(:);
b=b(:);

% lm gets initialized at the limit case in which a=b

lm = a;

tol = 1e-8;

%The cases in which b != a get distinguished and the updated according to
%the formula
idx = abs(b-a) > tol;

lm(idx) = (b(idx) - a(idx)) ./ log(b(idx)./a(idx));

end