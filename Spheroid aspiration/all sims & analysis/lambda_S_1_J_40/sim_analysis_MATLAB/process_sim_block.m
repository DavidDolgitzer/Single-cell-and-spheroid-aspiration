function [timeCommon, meanData, semData, lower_k, higher_k, mean_k, ...
          lower_tau, higher_tau, mean_tau, ...
          lower_xi_f, higher_xi_f, mean_xi_f, ...
          lower_xi_t, higher_xi_t, mean_xi_t, ...
          lower_L_p, higher_L_p, mean_L_p, ...
          lower_L_end, higher_L_end, mean_L_end, fitData] = ...
          process_sim_block(filePaths_asp, minLStart, asp_force, analysis_time_limit)

%PROCESS_SIM_BLOCK Process aspiration simulations with boundary correction.
%
%   Processes multiple simulation outputs, corrects boundary artifacts,
%   aligns valid corrected datasets to a common time grid, calculates
%   statistics, and fits a three-element viscoelastic model.
%
%   Inputs:
%   -------
%   filePaths_asp : cell array of strings
%       Paths to CSV files containing aspiration simulation outputs.
%       Each file must contain 'time' and 'L' columns.
%
%   minLStart : double, optional (default = 0)
%       Initial aspiration-length threshold to remove.
%
%   asp_force : double, required
%       Applied aspiration force in lambda-D units.
%
%   analysis_time_limit : double, optional (default = 20000)
%       Maximum analysis time in simulation-time units.
%
%   Outputs:
%   --------
%   timeCommon : array
%       Common time vector in seconds.
%
%   meanData : array
%       Mean aspiration length across valid simulations.
%
%   semData : array
%       Standard error of the mean aspiration length.
%
%   lower_k, higher_k, mean_k : double
%       Spring constant confidence bounds and fitted value.
%
%   lower_tau, higher_tau, mean_tau : double
%       Retardation-time confidence bounds and fitted value.
%
%   lower_xi_f, higher_xi_f, mean_xi_f : double
%       Fast-viscosity confidence bounds and fitted value.
%
%   lower_xi_t, higher_xi_t, mean_xi_t : double
%       Terminal-viscosity confidence bounds and fitted value.
%
%   lower_L_p, higher_L_p, mean_L_p : double
%       Plateau aspiration-length confidence bounds and fitted value.
%
%   lower_L_end, higher_L_end, mean_L_end : double
%       Aspiration length at the selected endpoint, reported using
%       mean plus or minus one SEM.
%
%   fitData : array
%       Fitted model evaluated at timeCommon.

% =========================================================================
% INPUT VALIDATION
% =========================================================================

if nargin < 2 || isempty(minLStart)
    minLStart = 0;
end

if nargin < 3 || isempty(asp_force)
    error('asp_force must be provided.');
end

if nargin < 4 || isempty(analysis_time_limit)
    analysis_time_limit = 20000;
end

if isempty(filePaths_asp)
    error('filePaths_asp must contain at least one file path.');
end

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% Time conversion: simulation time units to seconds
time_conversion = 45 / 7500;

% Geometric parameters for scaling
cell_area_nodes = 80;
num_cells = 500;
cell_radius_real = 5;

cell_area_real = pi * cell_radius_real^2;
total_area_real = cell_area_real * num_cells;
total_area_nodes = cell_area_nodes * num_cells;

spheroid_radius_real = sqrt(total_area_real / pi);
spheroid_radius_nodes = sqrt(total_area_nodes / pi);
length_scaling = spheroid_radius_real / spheroid_radius_nodes;

% Keep the external limit in simulation units and use seconds internally
analysis_time_limit_seconds = analysis_time_limit * time_conversion;

% Calculate effective force
baseline_force = 0.26;
force_effective = asp_force - baseline_force;

if force_effective <= 0
    error(['Effective aspiration force must be positive. ' ...
           'Received asp_force = %.4f and baseline force = %.4f.'], ...
          asp_force, baseline_force);
end

fprintf('Effective aspiration force: %.4f lambda-D\n', force_effective);

% =========================================================================
% BOUNDARY PARAMETERS
% =========================================================================

boundary_threshold = 130 * length_scaling - minLStart;

% =========================================================================
% DEFINE VISCOELASTIC MODEL
% =========================================================================

% Parameters:
%   params(1) = J
%   params(2) = phi
%   params(3) = tau
%   params(4) = t0

viscoelastic_model = @(params, t) ...
    force_effective * params(1) .* ...
    (1 - exp(-max(0, t - params(4)) ./ params(3))) + ...
    force_effective * params(2) .* max(0, t - params(4));

% =========================================================================
% LOAD AND PREPROCESS DATASETS
% =========================================================================

num_simulations = length(filePaths_asp);
processedDatasets = cell(num_simulations, 1);
valid_dataset = true(num_simulations, 1);

fprintf('Loading and preprocessing %d datasets...\n', num_simulations);

for i = 1:num_simulations
    data = readtable(filePaths_asp{i}, 'VariableNamingRule', 'preserve');

    if ~all(ismember({'time', 'L'}, data.Properties.VariableNames))
        error('File must contain time and L columns: %s', filePaths_asp{i});
    end

    % Find the first point exceeding the starting-length threshold
    first_valid = find(data.L > minLStart, 1, 'first');

    % Keep one point before aspiration begins when such a point exists
    if isempty(first_valid)
        start_index = 1;
    else
        start_index = max(first_valid - 1, 1);
    end

    data = data(start_index:end, :);

    % Reset time to zero and convert to seconds
    data.time = (data.time - data.time(1)) * time_conversion;

    % Restrict data to the selected analysis interval
    last_valid_time = find(data.time <= analysis_time_limit_seconds, 1, 'last');

    if isempty(last_valid_time)
        fprintf(['Warning: Dataset %d (%s) contains no data within the ' ...
                 'analysis interval and will be excluded.\n'], ...
                i, filePaths_asp{i});

        valid_dataset(i) = false;
        processedDatasets{i} = data;
        continue;
    end

    data = data(1:last_valid_time, :);

    % Scale length and apply the starting-length offset
    data.L = data.L * length_scaling - minLStart;

    processedDatasets{i} = data;
end

fprintf('Preprocessing complete.\n');

% =========================================================================
% DETECT BOUNDARY-HIT DATASETS
% =========================================================================

fprintf('Detecting boundary artifacts...\n');

hit_boundary = false(num_simulations, 1);
time_hit_boundary = NaN(num_simulations, 1);

for i = 1:num_simulations
    if ~valid_dataset(i)
        continue;
    end

    data = processedDatasets{i};
    boundary_index = find(data.L >= boundary_threshold, 1, 'first');

    if ~isempty(boundary_index)
        hit_boundary(i) = true;
        time_hit_boundary(i) = data.time(boundary_index);
    end
end

fprintf('Found %d valid datasets reaching the boundary threshold.\n', ...
        sum(hit_boundary & valid_dataset));

% =========================================================================
% CORRECT BOUNDARY ARTIFACTS
% =========================================================================

fprintf('Correcting boundary artifacts...\n');

correction_options = optimset( ...
    'Display', 'off', ...
    'MaxIter', 1000, ...
    'TolFun', 1e-8);

for i = 1:num_simulations
    if ~valid_dataset(i) || ~hit_boundary(i)
        continue;
    end

    data = processedDatasets{i};

    % Select points before boundary contact
    valid_for_fit = data.L < boundary_threshold & ...
                    data.time < time_hit_boundary(i);

    time_fit = data.time(valid_for_fit);
    length_fit = data.L(valid_for_fit);

    if length(time_fit) < 10
        fprintf(['Warning: Dataset %d (%s) has only %d usable ' ...
                 'pre-threshold points. At least 10 are required for ' ...
                 'boundary correction. The dataset will be excluded.\n'], ...
                i, filePaths_asp{i}, length(time_fit));

        valid_dataset(i) = false;
        continue;
    end

    % Initial estimates for boundary-correction fit
    J_guess = max(length_fit) / force_effective;

    midIdx = max(1, floor(length(time_fit) / 2));
    time_diff = time_fit(end) - time_fit(midIdx);

    if time_diff <= 0
        phi_guess = 1e-6;
    else
        phi_guess = max( ...
            (length_fit(end) - length_fit(midIdx)) / ...
            time_diff / force_effective, ...
            1e-6);
    end

    tau_guess = max(time_fit(midIdx), 1e-3);
    t0_guess = 0;

    initial_params = [J_guess, phi_guess, tau_guess, t0_guess];
    lower_bounds = [0, 0, 1e-3, -2];
    upper_bounds = [Inf, Inf, Inf, 5];

    try
        fitted_params = lsqcurvefit( ...
            viscoelastic_model, initial_params, time_fit, length_fit, ...
            lower_bounds, upper_bounds, correction_options);

        % Determine replacement region
        anchor_index = find(data.time <= time_hit_boundary(i), 1, 'last');
        anchor_length = data.L(anchor_index);

        replace_indices = data.time > time_hit_boundary(i);
        time_to_replace = data.time(replace_indices);

        % Validate the temporal sampling interval
        time_step = median(diff(data.time));

        if isempty(time_step) || ~isfinite(time_step) || time_step <= 0
            fprintf(['Warning: Dataset %d (%s) has no valid positive ' ...
                     'time step. Boundary extrapolation cannot be ' ...
                     'performed, and the dataset will be excluded.\n'], ...
                    i, filePaths_asp{i});

            valid_dataset(i) = false;
            continue;
        end

        % Create extrapolation times
        if max(data.time) < analysis_time_limit_seconds
            time_extrapolated = ((max(data.time) + time_step): ...
                                 time_step:analysis_time_limit_seconds)';
        else
            time_extrapolated = [];
        end

        % Apply boundary-correction model
        time_for_model = [time_to_replace; time_extrapolated];
        length_predicted = viscoelastic_model(fitted_params, time_for_model);

        % Prevent retraction and enforce monotonic aspiration
        if isempty(length_predicted)
            length_corrected = [];
        else
            length_corrected = max(length_predicted, anchor_length);
            length_corrected = cummax(length_corrected);
        end

        % Separate replacement and appended values
        num_replaced = length(time_to_replace);
        length_replaced = length_corrected(1:num_replaced);
        length_appended = length_corrected(num_replaced + 1:end);

        % Replace existing post-threshold values
        processedDatasets{i}.L(replace_indices) = length_replaced;

        % Append extrapolated values
        if ~isempty(time_extrapolated)
            extra_rows = processedDatasets{i}( ...
                repmat(height(processedDatasets{i}), ...
                       length(time_extrapolated), 1), :);

            extra_rows.time = time_extrapolated;
            extra_rows.L = length_appended;

            processedDatasets{i} = [processedDatasets{i}; extra_rows];
        end

    catch ME
        fprintf(['Warning: Boundary correction failed for dataset %d ' ...
                 '(%s): %s\n'], ...
                i, filePaths_asp{i}, ME.message);

        fprintf('Dataset %d will be excluded from all subsequent calculations.\n', i);
        valid_dataset(i) = false;
    end
end

fprintf('Boundary correction complete.\n');

% =========================================================================
% REPORT EXCLUDED DATASETS
% =========================================================================

excluded_indices = find(~valid_dataset);
valid_indices = find(valid_dataset);

if isempty(excluded_indices)
    fprintf('All %d datasets remain valid for ensemble analysis.\n', ...
            num_simulations);
else
    fprintf(['Excluded %d of %d datasets from ensemble analysis. ' ...
             'Excluded dataset indices: %s\n'], ...
            length(excluded_indices), num_simulations, ...
            mat2str(excluded_indices'));
end

if isempty(valid_indices)
    error(['No valid datasets remain after preprocessing and boundary ' ...
           'correction. Ensemble statistics cannot be calculated.']);
end

fprintf('Using %d valid datasets for ensemble analysis.\n', ...
        length(valid_indices));

% =========================================================================
% ALIGN DATA AND CALCULATE ENSEMBLE STATISTICS
% =========================================================================

fprintf('Computing ensemble statistics...\n');

num_time_points = 100;
timeCommon = linspace(0, analysis_time_limit_seconds, num_time_points);

num_valid_simulations = length(valid_indices);
aligned_lengths = NaN(num_valid_simulations, num_time_points);

for row_idx = 1:num_valid_simulations
    dataset_idx = valid_indices(row_idx);

    aligned_lengths(row_idx, :) = interp1( ...
        processedDatasets{dataset_idx}.time, ...
        processedDatasets{dataset_idx}.L, ...
        timeCommon, 'linear', NaN);
end

meanData = mean(aligned_lengths, 1, 'omitnan');
num_valid_at_time = sum(isfinite(aligned_lengths), 1);
semData = std(aligned_lengths, [], 1, 'omitnan') ./ sqrt(num_valid_at_time);

% =========================================================================
% CALCULATE LENGTH AT ANALYSIS ENDPOINT
% =========================================================================

target_time = analysis_time_limit_seconds;
% target_time = 120;

mean_L_end = interp1(timeCommon, meanData, target_time);
sem_L_end = interp1(timeCommon, semData, target_time);

lower_L_end = mean_L_end - sem_L_end;
higher_L_end = mean_L_end + sem_L_end;

fprintf('Aspiration length at %.2f s: %.2f +/- %.2f micrometres\n', ...
        target_time, mean_L_end, sem_L_end);

% =========================================================================
% PREPARE HIGH-RESOLUTION DATA FOR FITTING
% =========================================================================

fprintf('Computing high-resolution mean curve...\n');

num_fit_points = 1000;
time_fit = linspace(0, analysis_time_limit_seconds, num_fit_points);
aligned_fit = NaN(num_valid_simulations, num_fit_points);

for row_idx = 1:num_valid_simulations
    dataset_idx = valid_indices(row_idx);

    aligned_fit(row_idx, :) = interp1( ...
        processedDatasets{dataset_idx}.time, ...
        processedDatasets{dataset_idx}.L, ...
        time_fit, 'linear', NaN);
end

mean_fit = mean(aligned_fit, 1, 'omitnan');

time_fit_vec = time_fit(:);
length_fit_vec = mean_fit(:);

valid_mask = isfinite(time_fit_vec) & isfinite(length_fit_vec);
time_valid = time_fit_vec(valid_mask);
length_valid = length_fit_vec(valid_mask);

if length(time_valid) < 5
    error(['Insufficient valid ensemble points remain for the final ' ...
           'four-parameter model fit.']);
end

fprintf('Fitting %d data points...\n', length(time_valid));

% =========================================================================
% ESTIMATE INITIAL PARAMETERS
% =========================================================================

compliance_guess = max(length_valid) / force_effective;

mid_point = floor(length(length_valid) / 2);
time_difference = time_valid(end) - time_valid(mid_point);

if time_difference <= 0
    phi_guess = 0;
else
    phi_guess = (length_valid(end) - length_valid(mid_point)) / ...
                time_difference / force_effective;

    phi_guess = max(phi_guess, 0);
end

length_exponential = (length_valid / force_effective) - ...
                     phi_guess * time_valid;

length_exponential = length_exponential - min(length_exponential);

tau_index = find( ...
    length_exponential > 0.5 * max(length_exponential), ...
    1, 'first');

if isempty(tau_index)
    tau_guess = max(1e-3, time_valid(min(end, 2)));
else
    tau_guess = max(1e-6, time_valid(tau_index));
end

t0_guess = 0;

initial_params = [compliance_guess, phi_guess, tau_guess, t0_guess];
lower_bounds = [0, 0, 1e-3, -2];
upper_bounds = [Inf, Inf, Inf, 5];

fprintf('Initial parameter estimates:\n');
fprintf('  J   = %.4e micrometres/lambda-D\n', initial_params(1));
fprintf('  phi = %.4e micrometres/(lambda-D*s)\n', initial_params(2));
fprintf('  tau = %.4e s\n', initial_params(3));
fprintf('  t0  = %.4e s\n', initial_params(4));

% =========================================================================
% FIT MODEL TO ENSEMBLE MEAN
% =========================================================================

[fit_params, ~, residuals, ~, ~, ~, jacobian] = lsqcurvefit( ...
    viscoelastic_model, initial_params, time_valid, length_valid, ...
    lower_bounds, upper_bounds);

fprintf('\nFitted parameters:\n');
fprintf('  J   = %.4e micrometres/lambda-D\n', fit_params(1));
fprintf('  phi = %.4e micrometres/(lambda-D*s)\n', fit_params(2));
fprintf('  tau = %.4e s\n', fit_params(3));
fprintf('  t0  = %.4e s\n', fit_params(4));

% =========================================================================
% CALCULATE CONFIDENCE INTERVALS
% =========================================================================

conf_intervals = nlparci( ...
    fit_params, residuals, 'jacobian', jacobian);

num_observations = length(length_valid);
num_parameters = length(fit_params);
degrees_freedom = num_observations - num_parameters;

if degrees_freedom <= 0
    error('The final fit has insufficient degrees of freedom.');
end

residual_variance = sum(residuals.^2) / degrees_freedom;
param_covariance = residual_variance * ...
                   ((jacobian' * jacobian) \ eye(num_parameters));

% =========================================================================
% CALCULATE PHYSICAL PARAMETERS
% =========================================================================

% Spring constant: k = 1/J
mean_k = 1 / fit_params(1);
ci_k = [1 / conf_intervals(1, 2), 1 / conf_intervals(1, 1)];

lower_k = ci_k(1);
higher_k = ci_k(2);

% Retardation time
mean_tau = fit_params(3);
ci_tau = conf_intervals(3, :);

lower_tau = ci_tau(1);
higher_tau = ci_tau(2);

% Terminal viscosity: xi_t = 1/phi
mean_xi_t = 1 / fit_params(2);
ci_xi_t = [1 / conf_intervals(2, 2), 1 / conf_intervals(2, 1)];

lower_xi_t = ci_xi_t(1);
higher_xi_t = ci_xi_t(2);

% Fast viscosity: xi_f = tau/J
mean_xi_f = fit_params(3) / fit_params(1);

gradient_xi_f = [ ...
    -fit_params(3) / fit_params(1)^2, ...
    0, ...
    1 / fit_params(1), ...
    0];

se_xi_f = sqrt( ...
    gradient_xi_f * param_covariance * gradient_xi_f');

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
% PRINT SUMMARY
% =========================================================================

fprintf('\n========== SUMMARY OF RESULTS ==========\n');
fprintf('Valid datasets used: %d of %d\n', num_valid_simulations, num_simulations);
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

% =========================================================================
% EVALUATE FITTED MODEL
% =========================================================================

fitData = viscoelastic_model(fit_params, timeCommon(:));

end