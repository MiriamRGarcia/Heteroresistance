%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs for constant input to calculate extinction probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = Odes_Pext(tt, yy, mug, muk, AA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
%%% tt = Time,
%%% xx = CFUS,
%%% AA = Coefficient matrix,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
nr   = numel(mug);

% ----------------------------------------- %
% (1) Calculate total growth-death-net rates:
xx   = yy(1:nr);
dxdt = AA*xx;

xT   = sum(xx);

mugT = mug.'*xx/xT;
mukT = muk.'*xx/xT;

munT = mugT - mukT;

% ---------------------------------------- %
% (2) ODEs for the extinction probability:
zrho    = yy(nr + 1);

dzrhodt = - munT*zrho;

zz      = yy(nr + 2);

dzdt    = mukT*zrho;

% ---------------------------------------- %
% (3) Complete system of ODEs:
dydt    = [dxdt;
          dzrhodt;
          dzdt];

end