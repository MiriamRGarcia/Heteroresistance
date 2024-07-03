%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiBD_SSA: SSA (direct method) to simulate the multivariate 
%              BD heteroresistance model with
%              constant antimicrobial concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N, N_T] = multiBD_SSA(tsim, b, d, Xi, trans, N_0, N_TL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% tsim  = Sampling times to keep the state trajectories;
% b     = Array m_r x 1 with time dependent birth-rates for the m_r
%         different subpopulations at times tsim;
% d     = Array m_r x 1 with kill rates for the m_r subpopulations;
% Xi    = Matrix m_r x m_r with rates of change in AMR level;
% trans = Matrix m_r x {number of reactions} with the state transitions;
% N_0   = Array m_r x 1 of initial cell counts;
% N_TL  = Threshold value for total cell counts to stop simulation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Obtain problem sizes:
m_t = numel(tsim);
m_r = numel(N_0);

% ----------------------------------------------------------------------- %
% Initialice variables:

% Time:
t0         = tsim(1);                                                      % Initial time;
tf         = tsim(end);                                                    % Final time;
t          = t0;                                                           % Current simulation time;
t_intcount = 1;                                                            % Current simulation time interval;

% States:
N           = zeros(m_t, m_r);                                             % Matrix m_t x m_r of cell counts;
N(1, 1:m_r) = N_0.';                                                       % Save cell counts at t0;
N_T0        = sum(N_0);                                                    % Total cell count;
N_T         = sum(N, 2);                                                   % Initialise array m_t x 1 of total cell counts;

% Propensities:
Xi     = reshape(Xi.', [], 1);
Xi_aux = Xi.*reshape(repmat(N_0.', m_r, 1), [], 1);
Xi_aux(Xi_aux == 0) = [];

a  = [b.*N_0;
      d.*N_0;
      Xi_aux];
  
a_sum = sum(a);                                                            % Sum of propensities;


% ----------------------------------------------------------------------- %
% Main loop of SSA:

while t < tf 
    
    tau = -log(rand)/a_sum;                                                % Calculate time to the next reaction,
    t   = t + tau;                                                         % Actualice current time,
        
    % If the current time interval changes:
    if tsim(t_intcount + 1) < t
        t_intcount           = t_intcount + 1;                             % Actualice current interval counter;
        t                    = tsim(t_intcount);                           % Actualice time;
        N(t_intcount, 1:m_r) = N_0.';                                      % Actualice cell counts;
        N_T(t_intcount)      = N_T0;                                       % Actualice total cell counts;
        
        if m_t - 1 < t_intcount;return;end
    end
    
    % Calculate index of the next reaction:
    inext = find(cumsum(a) - rand*a_sum > 0, 1, 'first');                  
    
    % Actualice array of cell counts:
    N_0  = N_0 + trans(1:m_r, inext);
    
    % Actualice total cell counts:
    N_T0 = sum(N_0);                                 
    
    % Stop simulation if total counts reach the limit value:
    if N_T0 > N_TL        
        
        t_intcount               = t_intcount + 1;
        N(t_intcount:m_t, 1:m_r) = NaN;
        N_T(t_intcount:m_t)      = NaN;
        
        return
        
    % Stop simulation if the population goes extinct:  
    elseif N_T0 < 1 
        
        t_intcount               = t_intcount + 1;
        N(t_intcount:m_t, 1:m_r) = 0;
        N_T(t_intcount:m_t)      = 0;
        
        return
        
    end
    
    % Actualice propensities:
    Xi_aux = Xi.*reshape(repmat(N_0.', m_r, 1), [], 1);
    Xi_aux(Xi_aux == 0) = [];

    a  = [b.*N_0;
          d.*N_0;
          Xi_aux];
  
    a_sum = sum(a);                                                        % Sum of propensities;

end

% ----------------------------------------------------------------------- %
% If the time to next reaction exceds the final time:
N(t_intcount:m_t, 1:m_r) = repmat(N(t_intcount, 1:m_r), m_t - t_intcount + 1, 1);
N_T(t_intcount:m_t)      = repmat(N_T(t_intcount), m_t - t_intcount + 1);

end