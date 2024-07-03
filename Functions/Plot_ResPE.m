%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotPE: Plot results of model calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_ResPE(r, R, tsim, texp, Cexp, pars, logN_Tave_data, col, mks, ODEoptions, sim_name) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Problem sizes:
m_e       = numel(Cexp);
m_t       = numel(tsim);
m_texp    = numel(texp);

% Calculate total average counts:
[~, N_Tmod] = feval(sim_name, r, R, tsim, Cexp, pars, ODEoptions);
logN_Tmod   = log10(N_Tmod);

% Plot results:
figure

hold on

lgd  = cell(m_e, 1);

for iexp = 1:m_e
  
    plot(texp, logN_Tave_data(1:m_texp,iexp), mks{iexp}, 'Color', col(iexp,:), ...
        'MarkerSize', 10, 'MarkerFaceColor', col(iexp,:), 'MarkerEdgeColor',...
        'k', 'HandleVisibility', 'off')
    
    plot(tsim, logN_Tmod(1:m_t,iexp), 'Color', col(iexp, :), 'LineWidth', 1.5)
    
    lgd{iexp} = sprintf('$C=%0.2f$ (mg/L)', Cexp(iexp));

end

xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Total Bacterial Count ($\log_{10}N_T$)', 'Interpreter', 'Latex', 'FontSize', 10)

legend(lgd, 'Location', 'Best', 'Interpreter', 'Latex')

hold off

box off

end
