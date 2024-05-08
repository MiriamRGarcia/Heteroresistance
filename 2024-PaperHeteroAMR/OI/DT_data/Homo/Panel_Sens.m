%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel with sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close

% ----------------------------------------------------------------------- %
% Experimental setting:

% Parameter values (E.Coli):
mugS    = 0.63;%log(2)*3; %[h-1]
mugR    = 0.36;%0.96*mugS;
alph_g  = 2;

mukmaxS = 3.78;
bet     = 0.4;
alphk   = 3;

EC50k   = 1;
Hk      = 1;

xiSR    = 1e-6;
xiM     = 1e-2;
kxi     = log(1e2);

X0      = 1e6;
lamb_IC = 50;

par     = [mugS;mugR;alph_g;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lamb_IC];
    
np      = numel(par);

% Time discretisation:
t0       = 0; 
tf       = 48;       
ht       = 1e-2;             
tmod     = t0:ht:tf;        
texp     = [2 4 6 8 10 12 16 20 24 36 48];
tmod     = sort(unique([tmod texp])).'; 
texp_ind = find(ismember(tmod, texp));  
nt       = numel(tmod);    
ntexp    = numel(texp);

% Discretisation of the AMR level:
ra  = 0;
rb  = 1;
nr  = 50;
rr  = linspace(ra, rb, nr).';
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                              % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

% Parameters for input concentration:
MIC_S = EC50k*(mugS/(mukmaxS - mugS))^(1/Hk);
Cexp  = MIC_S*linspace(0, 8, 100);
Nexp  = numel(Cexp);

% Initial condition:
f0 = exp(-lamb_IC*rr);
f0 = f0/sum(f0);
x0 = X0*f0;

% ----------------------------------------------------------------------- %
% Calculate average and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

[xT, sens_xT] = Sens_MultiExp(tmod, rr, Cexp, par, x0, ODEoptions);

% ----------------------------------------------------------------------- %
% Plot sensitivities and normalised sensitivities if desired:

Plot_Sens(tmod, Cexp, par, xT, sens_xT)

%set(gcf, 'Position', [ -1456         116        1151         850])

% print('Panel_Fig2', '-depsc2')
% print('-r680','Panel_Fig2_HD', '-depsc2')