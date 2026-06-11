%%This script works as a pre-processing to the method application

%% Physical constants
q    = 1.602176634e-19;      % C
kB   = 1.380649e-23;         % J/K
T    = 300;                  % K
%In the project I will include the impact of hot electrons
VT   = kB*T/q;               % thermal voltage, V

c_lambda = 5/2;

%tau_E = [0.05, 0.1, 0.2, 0.5, 1.0]*1e-12
tau_E = 1e-12;

eps0  = 8.8541878128e-12;    % F/m
epsSi = 11.7*eps0;           % F/m

%% Device geometry
% L     = 1e-6;                % m
% Lc    = 0.15e-6;             % m
% sigma = 10e-9;               % m

L     = 3e-6;
Lc    = 0.45e-6; %thickness of high doped region
sigma = 30e-9; %this modulates the htan of the doping

%% Doping: 100 difference
ND_low  = 1e22;              % m^-3 = 1e16 cm^-3
ND_high = 1e24;              % m^-3 = 1e18 cm^-3

%% Electron mobility and diffusion coefficient
%also in the project mobility is constant
mun = 0.10;                  % m^2/(V s)
%NB: Einstein does not hold for differente temperatures
Dn  = mun*VT;                % m^2/s

%% Mesh-> uniform
N = 1000;                     % total number of nodes
x = linspace(0,L,N).';
h = x(2)-x(1);

Nint = N-2;                  % number of interior nodes

%% Smooth n+-n-n+ doping profile
% As given by the exercise
%First the slope: small sigma means more "step-like" doping, they are
%already shifted on the x axis
SL = 0.5*(1 - tanh((x - Lc)/sigma));
SR = 0.5*(1 + tanh((x - (L-Lc))/sigma));

ND = ND_low + (ND_high - ND_low).*(SL + SR);

%% Boundary conditions for density
n_left  = ND_high;
n_right = ND_high;

%% Boundary condition for electron thermal voltage theta
theta_left = VT;
theta_right = VT;

%% Bias continuation values
Va_list = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];  % V
%the device is symmetrical, thus i may apply negative voltages as well
%this set is used for the characteristic: the previous value is the initial
%guess of the next one

%Va_list = linspace(0,0.5,20);

% Va_list = 0;         % use this for equilibrium only

%% Useful scales
LD_low  = sqrt(epsSi*VT/(q*ND_low));
LD_high = sqrt(epsSi*VT/(q*ND_high));

%Voltage difference at equilibrium between high doped lateral zone and low 
% doped central zone
Vbi     = VT*log(ND_high/ND_low);

Jscale = (q*mun*VT*ND_high)/L;
fprintf('VT      = %.4g V\n', VT);
fprintf('Dn      = %.4g m^2/s\n', Dn);
fprintf('LD_low  = %.4g nm\n', LD_low*1e9);
fprintf('LD_high = %.4g nm\n', LD_high*1e9);
fprintf('Vbi     = %.4g V\n', Vbi);