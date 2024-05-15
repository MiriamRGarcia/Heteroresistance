%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OdesPext: Odes to calculate approximate extinction probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = OdesPext(t, z, b, d, AA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% t    = Time;
% z    = Auxiliary states;
% b    = Array of birth rates;
% d    = Array of death rates;
% AA   = Coefficient matrix;
%
% OUTPUT:
% dzdt = Time derivative of auxiliary states;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
nr = numel(b);

% ----------------------------------------- %
% Calculate total birth-death-net rates:
N    = z(1:nr);
dNdt = AA*N;

N_T  = sum(N);

b_T  = b.'*N/N_T;
d_T  = d.'*N/N_T;

g_T  = b_T - d_T;

% ---------------------------------------- %
% ODEs for the extinction probability:
zrho    = z(nr + 1);

dzrhodt = - g_T;

dzdt    = d_T*exp(zrho);

% ---------------------------------------- %
% Complete system of ODEs:
dzdt    = [dNdt;
          dzrhodt;
          dzdt];

end