function [timeCommon, meanData, semData, lower_k, higher_k, mean_k, ...
          lower_tau, higher_tau, mean_tau, ...
          lower_xi_f, higher_xi_f, mean_xi_f, ...
          lower_xi_t, higher_xi_t, mean_xi_t, ...
          lower_L_p, higher_L_p, mean_L_p, ...
          lower_L_end, higher_L_end, mean_L_end] = ...
          process_sim_block(filePaths_asp, clipping, asp_force)
%PROCESS_SIM_BLOCK Process aspiration simulation data and fit viscoelastic model
%
%   Processes multiple simulation outputs, aligns them to a common time grid,
%   calculates statistics, and fits a three-element viscoelastic model to
%   extract mechanical parameters (k, tau, xi_f, xi_t).
%
%   Inputs:
%   -------
%   filePaths_asp : cell array of strings
%       Paths to CSV files containing aspiration simulation outputs.
%       Each file must have 'time' and 'L' columns.
%   
%   clipping : double, optional (default = 0)
%       Initial length threshold to clip from the data (in simulation units).
%   
%   asp_force : double, required
%       Applied aspiration force (in λD units). Must be provided.
%
%   Outputs:
%   --------
%   timeCommon : array
%       Common time vector (in seconds) for aligned data
%   meanData : array
%       Mean aspiration length across simulations
%   semData : array
%       Standard error of the mean for aspiration length
%   lower_k, higher_k, mean_k : double
%       Spring constant confidence bounds and mean
%   lower_tau, higher_tau, mean_tau : double
%       Retardation time confidence bounds and mean
%   lower_xi_f, higher_xi_f, mean_xi_f : double
%       Fast viscosity confidence bounds and mean
%   lower_xi_t, higher_xi_t, mean_xi_t : double
%       Terminal viscosity confidence bounds and mean
%   lower_L_p, higher_L_p, mean_L_p : double
%       Plateau aspiration length confidence bounds and mean
%   lower_L_end, higher_L_end, mean_L_end : double
%       Aspiration length at 120s confidence bounds and mean
%   fitData : array
%       Fitted model evaluated at timeCommon points

% =========================================================================
% INPUT VALIDATION
% =========================================================================

% Set default clipping value if not provided
if nargin < 2 || isempty(clipping)
    clipping = 0;
end

% Aspiration force is required
if nargin < 3 || isempty(asp_force)
    error('asp_force must be provided.');
end

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% Time conversion: simulation time units → seconds
time_conversion = 45 / 7500;  % [s/simulation_unit]

% Geometric parameters for scaling
cell_area_nodes = 80;          % Area per cell in simulation [nodes]
num_cells = 500;               % Number of cells in spheroid
cell_radius_real = 5;          % Cell radius [μm]

% Calculate scaling factor from simulation to real units
cell_area_real = pi * cell_radius_real^2;
total_area_real = cell_area_real * num_cells;
total_area_nodes = cell_area_nodes * num_cells;

spheroid_radius_real = sqrt(total_area_real / pi);
spheroid_radius_nodes = sqrt(total_area_nodes / pi);

length_scaling = spheroid_radius_real / spheroid_radius_nodes;  % [μm/node]

% Calculate effective force (subtract baseline force)
baseline_force = 0.0;  % [λD] setting to zero since all we care about is tau
force_effective = asp_force - baseline_force;

fprintf('Effective aspiration force: %.4f λD\n', force_effective);

% =========================================================================
% LOAD AND PREPROCESS SIMULATION DATA
% =========================================================================

num_simulations = length(filePaths_asp);
sim_outputs = cell(num_simulations, 1);
max_time = Inf;  % Will be updated to shortest simulation duration

fprintf('Loading %d simulation files...\n', num_simulations);

for i = 1:num_simulations
    % Load CSV file
    data = readtable(filePaths_asp{i}, 'VariableNamingRule', 'preserve');
    
    % Find valid data range (after clipping threshold, before boundary)
    threshold = clipping;
    max_length = 159;  % Simulation boundary limit
    
    first_valid = find(data.L > threshold, 1, 'first');
    last_valid = find(data.L < max_length, 1, 'last');
    
    % Validate analysis interval
    if ~isempty(first_valid)
        start_idx = max(first_valid - 1, 1);
        
        if isempty(last_valid) || last_valid < start_idx
            error('No valid analysis interval in file: %s', filePaths_asp{i});
        end
        
        % Extract valid data range
        data = data(start_idx:last_valid, :);
        
        % Reset time to start at zero
        data.time = data.time - data.time(1);
    end
    
    % Scale length to physical units and apply clipping offset
    data.L = data.L * length_scaling - threshold;
    
    % Store processed data
    sim_outputs{i} = data;
    
    % Update maximum common time
    max_time = min(max_time, max(data.time));
end

fprintf('Data loaded. Common time range: 0 to %.2f simulation units\n', max_time);

% =========================================================================
% ALIGN DATA TO COMMON TIME GRID
% =========================================================================

% Create common time grid for all simulations
num_time_points = 100;
time_common_raw = linspace(0, max_time, num_time_points);

% Interpolate all simulations to common time grid
aligned_lengths = NaN(num_simulations, num_time_points);

for i = 1:num_simulations
    aligned_lengths(i, :) = interp1(...
        sim_outputs{i}.time, ...
        sim_outputs{i}.L, ...
        time_common_raw, ...
        'linear', ...
        NaN);
end

% =========================================================================
% CALCULATE STATISTICS
% =========================================================================

% Calculate mean and standard error across simulations
mean_length = mean(aligned_lengths, 1, 'omitnan');

num_valid = sum(isfinite(aligned_lengths), 1);
sem_length = std(aligned_lengths, [], 1, 'omitnan') ./ sqrt(num_valid);

% Convert time to seconds
timeCommon = time_common_raw * time_conversion;
meanData = mean_length;
semData = sem_length;

% =========================================================================
% CALCULATE ASPIRATION LENGTH AT 120 SECONDS
% =========================================================================

target_time = 120;  % [s]

mean_L_end = interp1(timeCommon, meanData, target_time);
sem_L_end = interp1(timeCommon, semData, target_time);

lower_L_end = mean_L_end - sem_L_end;
higher_L_end = mean_L_end + sem_L_end;

fprintf('Aspiration length at %.0f s: %.2f ± %.2f μm\n', ...
    target_time, mean_L_end, sem_L_end);

% =========================================================================
% PREPARE HIGH-RESOLUTION DATA FOR FITTING
% =========================================================================

% Create finer time grid for accurate fitting
num_fit_points = 1000;
time_fit_raw = linspace(0, max_time, num_fit_points);

% Interpolate to fine grid
aligned_fit = NaN(num_simulations, num_fit_points);

for i = 1:num_simulations
    aligned_fit(i, :) = interp1(...
        sim_outputs{i}.time, ...
        sim_outputs{i}.L, ...
        time_fit_raw, ...
        'linear', ...
        NaN);
end

% Calculate mean on fine grid
mean_fit = mean(aligned_fit, 1, 'omitnan');
time_fit = time_fit_raw * time_conversion;

% Extract valid (finite) data points for fitting
time_fit_vec = time_fit(:);
length_fit_vec = mean_fit(:);

valid_mask = isfinite(time_fit_vec) & isfinite(length_fit_vec);
time_valid = time_fit_vec(valid_mask);
length_valid = length_fit_vec(valid_mask);

fprintf('Fitting %d data points...\n', length(time_valid));

% =========================================================================
% DEFINE VISCOELASTIC MODEL
% =========================================================================

% Three-element viscoelastic model (Kelvin-Voigt + dashpot)
% L(t) = F * J * [1 - exp(-(t-t0)/tau)] + F * phi * (t-t0)
%
% Parameters:
%   params(1) = J     : instantaneous compliance [μm/λD]
%   params(2) = phi   : linear creep coefficient [μm/(λD·s)]
%   params(3) = tau   : retardation time [s]
%   params(4) = t0    : time offset [s]

viscoelastic_model = @(params, t) ...
    force_effective * params(1) * (1 - exp(-max(0, t - params(4)) / params(3))) + ...
    force_effective * params(2) * max(0, t - params(4));

% =========================================================================
% ESTIMATE INITIAL PARAMETERS
% =========================================================================

% Estimate instantaneous compliance from maximum length
compliance_guess = max(length_valid) / force_effective;

% Estimate linear creep from late-time slope
mid_point = floor(length(length_valid) / 2);
phi_guess = (...
    (length_valid(end) - length_valid(mid_point)) / ...
    (time_valid(end) - time_valid(mid_point))...
) / force_effective;
phi_guess = max(phi_guess, 0);  % Ensure non-negative

% Estimate retardation time from exponential rise
length_exponential = (length_valid / force_effective) - phi_guess * time_valid;
length_exponential = length_exponential - min(length_exponential);

tau_index = find(length_exponential > 0.5 * max(length_exponential), 1, 'first');
if isempty(tau_index)
    tau_guess = max(1e-3, time_valid(min(end, 2)));
else
    tau_guess = max(1e-6, time_valid(tau_index));
end

% Time offset estimate
t0_guess = 0;

% Compile initial guess
initial_params = [compliance_guess, phi_guess, tau_guess, t0_guess];

% Set parameter bounds
lower_bounds = [0,    0,    1e-3, -2];
upper_bounds = [Inf,  Inf,  Inf,   5];

fprintf('Initial parameter estimates:\n');
fprintf('  J   = %.4e μm/λD\n', initial_params(1));
fprintf('  phi = %.4e μm/(λD·s)\n', initial_params(2));
fprintf('  tau = %.4e s\n', initial_params(3));
fprintf('  t0  = %.4e s\n', initial_params(4));

% =========================================================================
% FIT MODEL TO DATA
% =========================================================================

[fit_params, ~, residuals, ~, ~, ~, jacobian] = lsqcurvefit(...
    viscoelastic_model, ...
    initial_params, ...
    time_valid, ...
    length_valid, ...
    lower_bounds, ...
    upper_bounds);

fprintf('\nFitted parameters:\n');
fprintf('  J   = %.4e μm/λD\n', fit_params(1));
fprintf('  phi = %.4e μm/(λD·s)\n', fit_params(2));
fprintf('  tau = %.4e s\n', fit_params(3));
fprintf('  t0  = %.4e s\n', fit_params(4));

% =========================================================================
% CALCULATE CONFIDENCE INTERVALS
% =========================================================================

% Standard 95% confidence intervals for fitted parameters
conf_intervals = nlparci(fit_params, residuals, 'jacobian', jacobian);

% Calculate degrees of freedom
num_observations = length(length_valid);
num_parameters = length(fit_params);
degrees_freedom = num_observations - num_parameters;

% Estimate residual variance
residual_variance = sum(residuals.^2) / degrees_freedom;

% Approximate parameter covariance matrix
% Using: cov(params) ≈ σ² * (J'J)^(-1)
param_covariance = residual_variance * ((jacobian' * jacobian) \ eye(num_parameters));

% =========================================================================
% CALCULATE PHYSICAL PARAMETERS WITH CONFIDENCE INTERVALS
% =========================================================================

% Spring constant: k = 1/J
mean_k = 1 / fit_params(1);
ci_k = [1/conf_intervals(1,2), 1/conf_intervals(1,1)];  % Reverse order
lower_k = ci_k(1);
higher_k = ci_k(2);

% Retardation time: tau (direct from fit)
mean_tau = fit_params(3);
ci_tau = conf_intervals(3, :);
lower_tau = ci_tau(1);
higher_tau = ci_tau(2);

% Terminal viscosity: xi_t = 1/phi
mean_xi_t = 1 / fit_params(2);
ci_xi_t = [1/conf_intervals(2,2), 1/conf_intervals(2,1)];  % Reverse order
lower_xi_t = ci_xi_t(1);
higher_xi_t = ci_xi_t(2);

% Fast viscosity: xi_f = tau/J (requires delta method)
mean_xi_f = fit_params(3) / fit_params(1);

% Gradient of xi_f with respect to [J, phi, tau, t0]
gradient_xi_f = [...
    -fit_params(3) / (fit_params(1)^2), ...  % ∂xi_f/∂J
    0, ...                                     % ∂xi_f/∂phi
    1 / fit_params(1), ...                    % ∂xi_f/∂tau
    0];                                        % ∂xi_f/∂t0

% Delta-method standard error: SE = sqrt(∇f' * Σ * ∇f)
se_xi_f = sqrt(gradient_xi_f * param_covariance * gradient_xi_f');

% 95% confidence interval
t_critical = tinv(0.975, degrees_freedom);
ci_xi_f = mean_xi_f + [-1, 1] * t_critical * se_xi_f;
lower_xi_f = ci_xi_f(1);
higher_xi_f = ci_xi_f(2);

% Plateau aspiration length: L_p = F_eff * J
mean_L_p = force_effective * fit_params(1);
ci_L_p = force_effective * conf_intervals(1, :);
lower_L_p = ci_L_p(1);
higher_L_p = ci_L_p(2);

% =========================================================================
% PRINT RESULTS SUMMARY
% =========================================================================

fprintf('\n=== FITTED MECHANICAL PARAMETERS ===\n');
fprintf('Spring constant (k):\n');
fprintf('  Mean: %.4e, 95%% CI: [%.4e, %.4e]\n', mean_k, lower_k, higher_k);
fprintf('Retardation time (tau):\n');
fprintf('  Mean: %.4e s, 95%% CI: [%.4e, %.4e] s\n', mean_tau, lower_tau, higher_tau);
fprintf('Fast viscosity (xi_f):\n');
fprintf('  Mean: %.4e, 95%% CI: [%.4e, %.4e]\n', mean_xi_f, lower_xi_f, higher_xi_f);
fprintf('Terminal viscosity (xi_t):\n');
fprintf('  Mean: %.4e, 95%% CI: [%.4e, %.4e]\n', mean_xi_t, lower_xi_t, higher_xi_t);
fprintf('Plateau length (L_p):\n');
fprintf('  Mean: %.4f μm, 95%% CI: [%.4f, %.4f] μm\n', mean_L_p, lower_L_p, higher_L_p);
fprintf('\n=====================================\n');
fprintf('\n=====================================\n');

end