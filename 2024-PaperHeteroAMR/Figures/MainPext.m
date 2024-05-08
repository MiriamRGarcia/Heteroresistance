%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for plot extinction probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

addpath('Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Define user settings:

% ------------------------------------------ %
% Parameter values:
bS        = 0.63;
bR        = 0.36;
alpha_b   = 2;

d_maxS    = 3.78;
beta_d    = 0.4;
alpha_d   = 3;

EC_50d    = 1;
H_d       = 1;

xi_SR     = 1e-6;
k_xi      = log(1e2);

N_T0      = 1e6;
lambda_T0 = 50;

par       = [bS;bR;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;xi_SR;k_xi;N_T0;lambda_T0];
    
np        = numel(par);

% Time discretisation:
t0   = 0;
tf   = 48;
ht   = 1e-3;
tmod = (t0:ht:tf).';
nt   = numel(tmod);

% Discretisation of AMR level:
nr   = 50;
ra   = 0;
rb   = 1;
hr   = (rb - ra)/(nr - 1);
rr   = (ra:hr:rb).';
RR   = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                      
RR   = RR - triu(RR) + tril(RR).';

% ------------------------------------------ %
% Antimicrobial concentration:
MIC_S = EC_50d*(bS/(d_maxS - bS))^(1/H_d);
C     = MIC_S*linspace(0, 64, 100);
nC    = numel(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Calculate extinction probability:

% Initial condition for mean cell numbers:
f0  = exp(-lambda_T0*rr);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice total population size:
NT   = zeros(nt, nC);

% Initialise extinction probability:
Pext = zeros(nt, nC);

% Initialise coefficient matrix:
Xi   = xi_SR*exp(k_xi*(1 - RR));
Xi   = Xi - diag(diag(Xi));

AA_aux = Xi.' - diag(sum(Xi, 2));

% Birth and maximal kill rates:
b     = bS*bR./(bR + rr.^alpha_b*(bS - bR));
d_max = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);

% Initial condition for ODEs:
z0 = [N0;
      0;
      0];

for iC = 1:nC
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = C(iC);
    
    % Calculate kill rate at current time:
    d_max  = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d      = d_max*HC;

    AA     = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Call to ODEs with extinction probability:
    [~, zout] = ode15s(@(t,s) OdesPext(t, s, b, d, AA), tmod, z0, ODEoptions);
    zz        = zout(1:nt, nr + 2);
    
    % Extinction probability:
    Pext(1:nt, iC) = (zz./(1 + zz)).^N_T0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Plot results (contour plot):

fig = figure;

set(gcf,'color','w');

fig.Position = [582   454   640   481];

Cticks = linspace(C(1), C(nC), 4);

% Call to contour plot:
contourf(tmod, C, Pext.', 50, 'edgecolor', 'none')

% Colorbar of subplot:
CB   = colorbar;
%lCB  = get(CB,'Limits');
%tCB  = linspace(lCB(1), lCB(2), 4);
    
CB.FontSize = 11;
    
%set(CB,'Ticks',tCB,'TickLabelInterpreter', 'Latex')
     
%TLCB = arrayfun(@(x) sprintf('%.2g',x), tCB, 'un', 0);
      
%set(CB,'TickLabels',TLCB)
set(CB, 'TickLabelInterpreter', 'Latex')
 
% Figure settings:
xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('$C$ (mg/L)', 'Interpreter','Latex','FontSize', 12);

xlim([tmod(1) tmod(end)])
ylim([C(1) C(end)])

% Axis properties:
set(gca,'XTick',[0 10 20 30 40 48], 'YTick', Cticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)

% Set ticks on y axis:
ax = gca;
ax.YAxis.TickLabelFormat = '%.2f';
    
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Print figure if desired:
print('-r720','PanelPextHD', '-dpng')


rmpath('Functions')