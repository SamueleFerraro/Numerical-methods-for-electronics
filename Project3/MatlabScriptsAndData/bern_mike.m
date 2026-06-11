function [bp,bn] = bern_mike(x)

% Vector implementation of the Bernoulli function.
% Returns bp = B(x) and bn = B(-x) for an input vector x.
%
% Observations:
%   for x >  709, exp( x) = inf ==> b returns  0
%   for x < -745, exp(-x) = 0   ==> b returns -x
%
% Properties:
%   bn = exp(x)*bp;
%   bn = x + bp;
%
%   Points of strength: 
%   - ensures bp(x) = bn(-x) for every x
%   - precision in on the order of epsM
%   - xlim = 2^-2 ensures maximum precision for every x 
%     at least for Taylor to the 10° order (trial-and-error)
% 

bp = zeros(size(x));
bn = bp;
xlim = 2^-2;

% |x| big & x positive
ibig         =  find(x > xlim);
if ~isempty(ibig)
  bp(ibig)   =  x(ibig)./(exp(x(ibig))-1);
  bn(ibig)   =  x(ibig) + bp(ibig);
end  

% |x| big & x negative
ibig         =  find(x < xlim);
if ~isempty(ibig)
  bn(ibig)   = -x(ibig)./(exp(-x(ibig))-1);
  bp(ibig)   =  bn(ibig) - x(ibig);
end

% |x| small
ismall       = find(abs(x) <= xlim);
if ~isempty(ismall)
  bp(ismall) = x(ismall).^10/47900160 - x(ismall).^8/1209600 + ...
               x(ismall).^6 /30240    - x(ismall).^4/720     + ...
               x(ismall).^2 /12       - x(ismall)   /2       + 1;

  bn(ismall) = x(ismall).^10/47900160 - x(ismall).^8/1209600 + ...
               x(ismall).^6 /30240    - x(ismall).^4/720     + ...
               x(ismall).^2 /12       + x(ismall)   /2       + 1;
end