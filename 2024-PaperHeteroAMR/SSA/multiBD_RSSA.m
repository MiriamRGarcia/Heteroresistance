%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiBD_RSSA: Rejection-based SSA to simulate the
%               multivariate BD heteroresistance model with
%               constant antimicrobial concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N, N_T] = multiBD_RSSA(tsim, b, d, Xi, trans, N_0, N_TL)
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
t_intcount = 1;                                                            % Current simulation time-interval;

% States:
N    = zeros(m_t, m_r);                                                    % Matrix m_t x m_r of cell counts;
N(1, 1:m_r) = N_0.';
N_T0 = sum(N_0);                                                           % Total cell count;
N_T  = sum(N, 2);                                                             % Array m_t x 1 of total cell counts;                                       

% Bounds on the current state value:
N_low = 0.85*N_0;
N_up  = 1.15*N_0;

% Bounds on propensities:
Xi    = reshape(Xi.', [], 1);

Xi_aux_low = Xi.*reshape(repmat(N_low.', m_r, 1), [], 1);
Xi_aux_low(Xi_aux_low == 0) = [];

a_low = [b.*N_low;
         d.*N_low;
         Xi_aux_low];
     
Xi_aux_up = Xi.*reshape(repmat(N_up.', m_r, 1), [], 1);
Xi_aux_up(Xi_aux_up == 0) = [];

a_up  = [b.*N_up;
         d.*N_up;
         Xi_aux_up];

a0_up = sum(a_up);

% ----------------------------------------------------------------------- %
% Main loop of tRSSA: 
while t < tf
    
    % Actualice current time:
    tau = - log(rand)/a0_up;
    t   = t + tau;
    
    % If the current time interval changes:
    if tsim(t_intcount + 1) < t
        t_intcount           = t_intcount + 1;                             % Actualice current interval counter;
        t                    = tsim(t_intcount);                           % Actualice time;
        N(t_intcount, 1:m_r) = N_0.';                                      % Actualice cell counts;
        N_T(t_intcount)      = N_T0;                                       % Actualice total cell counts;
        
        if m_t - 1 < t_intcount;return;end
    end
    
    % Acceptance step:  
    un1    = rand;
    un2    = rand;
    
    min_ind = find(cumsum(a_up) - un1*a0_up > 0, 1, 'first');
    
    if un2*a_up(min_ind) > a_low(min_ind)                                  % First acceptance condition;
        
        Xi_aux = Xi.*reshape(repmat(N_0.', m_r, 1), [], 1);
        Xi_aux(Xi_aux == 0) = [];
        
        a = [b.*N_0;
             d.*N_0;
             Xi_aux];
         
        a = a(min_ind);
        
        accept = 1;
        
        if un2*a_up(min_ind) > a
            accept = 0;
        end
        
    else                                                                   % Second acceptance condition;
        accept = 1;
    end
    
    if accept
   
        % Update current state:
        N_0 = N_0 + trans(1:m_r, min_ind);
        
        % Update total cell count:
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
        
        % Actualice bounds if necessary:
        if numel(find(N_0 - N_low < 0)) + numel(find(N_up - N_0 < 0)) > 0
            N_low = 0.85*N_0;
            N_up  = 1.15*N_0;
            
            Xi_aux_low = Xi.*reshape(repmat(N_low.', m_r, 1), [], 1);
            Xi_aux_low(Xi_aux_low == 0) = [];
            
            Xi_aux_up = Xi.*reshape(repmat(N_up.', m_r, 1), [], 1);
            Xi_aux_up(Xi_aux_up == 0) = [];
            
            a_low = [b.*N_low;
                     d.*N_low;
                     Xi_aux_low];
         
            a_up  = [b.*N_up;
                     d.*N_up;
                     Xi_aux_up];
                 
            a0_up = sum(a_up);
        end
    end
end

% ----------------------------------------------------------------------- %
% If the time to next reaction exceds the final time:
N(t_intcount:m_t, 1:m_r) = repmat(N(t_intcount, 1:m_r), m_t - t_intcount + 1, 1);
N_T(t_intcount:m_t)      = repmat(N_T(t_intcount), m_t - t_intcount + 1);


end