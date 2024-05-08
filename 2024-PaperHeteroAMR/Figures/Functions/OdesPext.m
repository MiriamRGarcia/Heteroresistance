%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODEs to calculate extinction probability (constant concentration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = OdesPext(tt, zz, b, d, AA)

% Problem sizes:
nr = numel(b);

% ----------------------------------------- %
% Calculate total birth-death-net rates:
N    = zz(1:nr);
dNdt = AA*N;

NT   = sum(N);

bT   = b.'*N/NT;
dT   = d.'*N/NT;

gT   = bT - dT;

% ---------------------------------------- %
% ODEs for the extinction probability:
zrho    = zz(nr + 1);

dzrhodt = - gT;

dzdt    = dT*exp(zrho);

% ---------------------------------------- %
% Complete system of ODEs:
dzdt    = [dNdt;
          dzrhodt;
          dzdt];

end