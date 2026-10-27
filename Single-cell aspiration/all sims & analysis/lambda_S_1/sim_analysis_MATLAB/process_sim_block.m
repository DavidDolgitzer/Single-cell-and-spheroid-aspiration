function [timeCommon, meanData, semData, lower_k, higher_k, mean_k, ...
          lower_tau, higher_tau, mean_tau, lower_xi_f, higher_xi_f, mean_xi_f, ...
          lower_L_p, higher_L_p, mean_L_p] = ...
          process_sim_block(filePaths_asp, asp_force)
%PROCESS_SIM_BLOCK Process aspiration simulation data and fit Kelvin-Voigt model
%
%   Processes multiple simulation outputs, aligns them to a common time grid,
%   calculates statistics, and fits a two-element Kelvin-Voigt viscoelastic 
%   model to extract mechanical parameters (k, tau, xi_f).
%
%   Inputs:
%   -------
%   filePaths_asp : cell array of strings
%       Paths to CSV files containing aspiration simulation outputs.
%       Each file must have 'time' and 'L' columns.
%   
%   asp_force : double, required
%       Applied aspiration force (in appropriate units).
%
%   Outputs:
%   --------
%   timeCommon : array
%       Common time vector (in seconds) for aligned data
%   meanData : array
%       Mean aspiration length across simulations [μm]
%   semData : array
%       Standard error of the mean for aspiration length [μm]
%   lower_k, higher_k, mean_k : double
%       Spring constant confidence bounds and mean
%   lower_tau, higher_tau, mean_tau : double
%       Retardation time confidence bounds and mean [s]
%   lower_xi_f, higher_xi_f, mean_xi_f : double
%       Friction coefficient confidence bounds and mean
%   lower_L_p, higher_L_p, mean_L_p : double
%       Plateau aspiration length confidence bounds and mean [μm]

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% Time conversion: simulation time units → seconds
time_conversion = 45 / 7500;  % 1 time step = 0.006 seconds

% Geometric parameters for length scaling
cell_area_nodes = 800;                        % Cell area in simulation [nodes²]
cell_radius_nodes = sqrt(cell_area_nodes/pi); % Cell radius [nodes]
cell_radius_real = 5;                         % Target cell radius [μm]

% Calculate scaling factor: simulation units → microns
length_scaling = cell_radius_real / cell_radius_nodes;  % [μm/node]

fprintf('Processing %d simulation files...\n', length(filePaths_asp));
fprintf('Length scaling factor: %.4f μm/node\n', length_scaling);

% =========================================================================
% LOAD AND PREPROCESS SIMULATION DATA
% =========================================================================

num_simulations = length(filePaths_asp);
sim_outputs = cell(num_simulations, 1);
max_time = Inf;  % Will be updated to shortest simulation duration

for i = 1:num_simulations
    % Load CSV file
    data = readtable(filePaths_asp{i}, 'VariableNamingRule', 'preserve');
    
    % Find first point where aspiration begins (L > 0)
    first_valid = find(data.L > 0, 1, 'first');
    
    if ~isempty(first_valid)
        % Keep one point before aspiration starts for smooth curve
        start_idx = max(first_valid - 1, 1);
        data = data(start_idx:end, :);
        
        % Reset time to start at zero
        data.time = data.time - data.time(1);
    end
    
    % Scale length to physical units (microns)
    data.L = data.L * length_scaling;
    
    % Store processed data
    sim_outputs{i} = data;
    
    % Update maximum common time across all simulations
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

% Calculate mean across all simulations
mean_length = mean(aligned_lengths, 1, 'omitnan');

% Calculate standard error of the mean (SEM)
% SEM accounts for the number of valid (non-NaN) observations at each time point
num_valid = sum(isfinite(aligned_lengths), 1);
sem_length = std(aligned_lengths, [], 1, 'omitnan') ./ sqrt(num_valid);

% Convert time to seconds for output
timeCommon = time_common_raw * time_conversion;
meanData = mean_length;
semData = sem_length;

fprintf('Statistics calculated over %d time points\n', num_time_points);

% =========================================================================
% PREPARE HIGH-RESOLUTION DATA FOR FITTING
% =========================================================================

% Create finer time grid for accurate curve fitting
num_fit_points = 10000;
time_fit_raw = linspace(0, max_time, num_fit_points);

% Interpolate all simulations to fine grid
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
length_fit_vec = mean_fit(:);  % [μm]

valid_mask = isfinite(time_fit_vec) & isfinite(length_fit_vec);
time_valid = time_fit_vec(valid_mask);
length_valid = length_fit_vec(valid_mask);

fprintf('Fitting %d high-resolution data points...\n', length(time_valid));

% =========================================================================
% DEFINE KELVIN-VOIGT MODEL
% =========================================================================

% Two-element Kelvin-Voigt model:
% L(t) = F * A * [1 - exp(-t/tau)]
%
% Where:
%   L(t) = aspiration length at time t [μm]
%   F    = aspiration force (constant)
%   A    = compliance (1/k) [μm/force_unit]
%   tau  = retardation time [s]
%
% Parameters:
%   params(1) = A   : compliance (inverse of spring constant)
%   params(2) = tau : retardation time

kelvin_voigt_model = @(params, t) ...
    asp_force * params(1) * (1 - exp(-t / params(2)));

% =========================================================================
% ESTIMATE INITIAL PARAMETERS
% =========================================================================

% Compliance estimate: plateau value divided by force
% At long times, L → F*A, so A ≈ L_max / F
compliance_guess = max(length_valid) / asp_force;

% Retardation time estimate: time to reach ~63.2% of maximum
% In exponential approach, 63.2% corresponds to one time constant
tau_threshold = 0.632 * max(length_valid);
tau_index = find(length_valid > tau_threshold, 1, 'first');

if isempty(tau_index)
    % If threshold not reached, use median time as fallback
    tau_guess = median(time_valid);
else
    tau_guess = time_valid(tau_index);
end

% Compile initial parameter guess
initial_params = [compliance_guess, tau_guess];

% Set parameter bounds (must be positive)
lower_bounds = [1e-12, 1e-12];
upper_bounds = [Inf, Inf];

fprintf('Initial parameter estimates:\n');
fprintf('  Compliance A = %.4e μm/force_unit\n', initial_params(1));
fprintf('  Tau          = %.4e s\n', initial_params(2));

% =========================================================================
% FIT MODEL TO DATA
% =========================================================================

[fit_params, ~, residuals, ~, ~, ~, jacobian] = lsqcurvefit(...
    kelvin_voigt_model, ...
    initial_params, ...
    time_valid, ...
    length_valid, ...
    lower_bounds, ...
    upper_bounds);

fprintf('\nFitted parameters:\n');
fprintf('  Compliance A = %.4e μm/force_unit\n', fit_params(1));
fprintf('  Tau          = %.4e s\n', fit_params(2));

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

% Spring constant: k = 1/A
mean_k = 1 / fit_params(1);
ci_k = [1/conf_intervals(1,2), 1/conf_intervals(1,1)];  % Reverse order for 1/A
lower_k = ci_k(1);
higher_k = ci_k(2);

% Retardation time: tau (direct from fit)
mean_tau = fit_params(2);
ci_tau = conf_intervals(2, :);
lower_tau = ci_tau(1);
higher_tau = ci_tau(2);

% Friction coefficient: xi_f = tau/A (requires delta method)
mean_xi_f = fit_params(2) / fit_params(1);

% Gradient of xi_f with respect to [A, tau]
gradient_xi_f = [...
    -fit_params(2) / (fit_params(1)^2), ...  % ∂xi_f/∂A
    1 / fit_params(1)];                      % ∂xi_f/∂tau

% Delta-method standard error: SE = sqrt(∇f' * Σ * ∇f)
se_xi_f = sqrt(gradient_xi_f * param_covariance * gradient_xi_f');

% 95% confidence interval using t-distribution
t_critical = tinv(0.975, degrees_freedom);
ci_xi_f = mean_xi_f + [-1, 1] * t_critical * se_xi_f;
lower_xi_f = ci_xi_f(1);
higher_xi_f = ci_xi_f(2);

% Plateau aspiration length: L_p = F * A
mean_L_p = asp_force * fit_params(1);
ci_L_p = asp_force * conf_intervals(1, :);
lower_L_p = ci_L_p(1);
higher_L_p = ci_L_p(2);

% =========================================================================
% PRINT RESULTS SUMMARY
% =========================================================================

fprintf('\n========== FITTED VISCOELASTIC PARAMETERS (95%% CI) ==========\n');
fprintf('Spring constant (k):\n');
fprintf('  %.4g  [%.4g, %.4g]\n', mean_k, lower_k, higher_k);
fprintf('\nRetardation time (tau):\n');
fprintf('  %.4g s  [%.4g, %.4g] s\n', mean_tau, lower_tau, higher_tau);
fprintf('\nFriction coefficient (xi_f):\n');
fprintf('  %.4g  [%.4g, %.4g]\n', mean_xi_f, lower_xi_f, higher_xi_f);
fprintf('\nPlateau aspiration length (L_p):\n');
fprintf('  %.4g μm  [%.4g, %.4g] μm\n', mean_L_p, lower_L_p, higher_L_p);
fprintf('\n=====================================\n');
fprintf('\n=====================================\n');

end