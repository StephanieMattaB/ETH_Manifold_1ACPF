%% RUN_IEEE13_NONLINEAR_AC
% Nonlinear AC validation for the modified three-phase IEEE 13-node feeder.
%
% Required existing file:
%   ieee13.m
%
% Required new file:
%   solve_acpf_3ph_rect.m
%
% The script solves:
%   1) nominal nonlinear AC power flow;
%   2) 1.8x stressed nonlinear AC power flow;
%   3) optionally, stressed + game-theoretic action.
%
% IMPORTANT CONVENTIONS
% ---------------------
% ieee13.m stores p and q as POSITIVE NET CONSUMPTION.
% solve_acpf_3ph_rect.m expects POSITIVE NETWORK INJECTION.
% Therefore:
%
%       Sspec = -Sconsumption.
%
% The fixed capacitor banks are separated before applying the 1.8 load
% multiplier. Otherwise, multiplying the returned q vector directly would
% incorrectly scale the capacitor injections too.

%clearvars;
%close all;
%clc;

%% Configuration

stress_factor = 1.3;
v_min = 0.95;
v_max = 1.05;

solver_opts = struct();
solver_opts.tol = 1e-10;
solver_opts.max_iter = 50;
solver_opts.max_line_search = 15;
solver_opts.verbose = true;

% Optional game result.
%
% Preferred MAT-file variables:
%   delta_s_ns : 72x1 real vector [delta_p_ns; delta_q_ns]
%
% Alternative MAT-file variables:
%   u_star     : physical game vector
%   game_cols  : column indices locating u_star inside the 72x1 vector
%
% game_convention:
%   'consumption' -> positive delta p/q means increased consumption
%   'injection'   -> positive delta p/q means increased network injection
game_file = '';                % Example: 'selected_gne_for_ac.mat'
game_convention = 'consumption';

%% Load the existing 1ACPF IEEE-13 representation

[Y, v_data, theta_data, p_net, q_net, n_bus] = ieee13();

Sbase = 5e6;
n_phase = 3;
N = n_phase*n_bus;

slack_idx = (1:3).';
ns_idx = (4:N).';

if ~isequal(size(Y), [N, N])
    error('Unexpected Y dimension. Expected %d-by-%d.', N, N);
end

Vslack = v_data(slack_idx) .* exp(1i*theta_data(slack_idx));

% ieee13.m returns positive net consumption:
%   S_net_consumption = S_load - S_capacitor.
S_net_nominal = p_net + 1i*q_net;
S_net_nominal(slack_idx) = 0;

%% Separate fixed capacitor injections from nominal demand

% Terminal index helper: phases are ordered [a,b,c] at each bus.
terminal = @(bus, phase) 3*(bus - 1) + phase;

S_cap = zeros(N, 1);

% ieee13.m:
%   bus 9:  200 kvar on phases a,b,c
%   bus 12: 100 kvar on phase c
S_cap(terminal(9, 1:3)) = 1i*(200e3/Sbase);
S_cap(terminal(12, 3)) = 1i*(100e3/Sbase);

% Since S_net_nominal = S_load_nominal - S_cap:
S_load_nominal = S_net_nominal + S_cap;

reconstruction_error = norm( ...
    (S_load_nominal - S_cap) - S_net_nominal, inf);

if reconstruction_error > 1e-14
    error('Load/capacitor decomposition failed.');
end

%% Case 1: nominal nonlinear AC power flow

% Convert positive consumption to positive network injection.
Sspec_nominal = -S_net_nominal;

% Flat balanced initialization using the fixed slack phasors.
V0 = repmat(Vslack, n_bus, 1);

fprintf('\n============================================================\n');
fprintf('CASE 1: NOMINAL NONLINEAR AC POWER FLOW\n');
fprintf('============================================================\n');

nominal = solve_acpf_3ph_rect( ...
    Y, slack_idx, Vslack, Sspec_nominal, V0, solver_opts);

assert_converged(nominal, 'nominal');

%% Compare nominal nonlinear result with published IEEE feeder voltages

v_testfeeder = [ ...
    1.0625  1.0500  1.0687; ...
    1.0210  1.0420  1.0174; ...
    1.0180  1.0401  1.0148; ...
    0.9940  1.0218  0.9960; ...
    NaN     1.0329  1.0155; ...
    NaN     1.0311  1.0134; ...
    0.9900  1.0529  0.9778; ...
    0.9900  1.0529  0.9777; ...
    0.9835  1.0553  0.9758; ...
    0.9900  1.0529  0.9778; ...
    0.9881  NaN     0.9758; ...
    NaN     NaN     0.9738; ...
    0.9825  NaN     NaN ...
];

v_test_vector = reshape(v_testfeeder.', N, 1);
benchmark_mask = ~isnan(v_test_vector);

benchmark_error = nominal.v(benchmark_mask) ...
    - v_test_vector(benchmark_mask);

nominal.benchmark_max_error = max(abs(benchmark_error));
nominal.benchmark_rms_error = sqrt(mean(benchmark_error.^2));

fprintf('\nNominal benchmark comparison:\n');
fprintf('  max |V_AC - V_benchmark| = %.6e p.u.\n', ...
    nominal.benchmark_max_error);
fprintf('  RMS |V_AC - V_benchmark| = %.6e p.u.\n', ...
    nominal.benchmark_rms_error);

%% Case 2: 1.8x stressed nonlinear AC power flow

% Scale physical demand only; keep capacitor banks fixed.
S_net_stressed = stress_factor*S_load_nominal - S_cap;
Sspec_stressed = -S_net_stressed;

fprintf('\n============================================================\n');
fprintf('CASE 2: %.2fx STRESSED NONLINEAR AC POWER FLOW\n', ...
    stress_factor);
fprintf('============================================================\n');

stressed = solve_acpf_3ph_rect( ...
    Y, slack_idx, Vslack, Sspec_stressed, nominal.V, solver_opts);

assert_converged(stressed, 'stressed');

%% Case 3: optional stressed + game action

game = [];
has_game_case = ~isempty(game_file);

if has_game_case
    if exist(game_file, 'file') ~= 2
        error('Game file not found: %s', game_file);
    end

    game_data = load(game_file);
    n_s_ns = 2*numel(ns_idx);   % [delta_p_ns; delta_q_ns] = 72x1

    if isfield(game_data, 'delta_s_ns')
        delta_s_ns = game_data.delta_s_ns(:);

        if numel(delta_s_ns) ~= n_s_ns
            error('delta_s_ns must have %d entries.', n_s_ns);
        end

    elseif isfield(game_data, 'u_star') && isfield(game_data, 'game_cols')
        u_star = game_data.u_star(:);
        game_cols = game_data.game_cols(:);

        if numel(u_star) ~= numel(game_cols)
            error('u_star and game_cols must have equal lengths.');
        end
        if any(game_cols < 1) || any(game_cols > n_s_ns)
            error('game_cols contains an invalid index.');
        end
        if numel(unique(game_cols)) ~= numel(game_cols)
            error('game_cols contains duplicate indices.');
        end

        delta_s_ns = zeros(n_s_ns, 1);
        delta_s_ns(game_cols) = u_star;

    else
        error([ ...
            'The MAT-file must contain either delta_s_ns, or both ', ...
            'u_star and game_cols.' ...
        ]);
    end

    n_ns = numel(ns_idx);
    delta_p_ns = delta_s_ns(1:n_ns);
    delta_q_ns = delta_s_ns(n_ns + 1:end);

    delta_S_full = zeros(N, 1);
    delta_S_full(ns_idx) = delta_p_ns + 1i*delta_q_ns;

    switch lower(game_convention)
        case 'consumption'
            % Positive game action means more consumption, hence less injection.
            delta_S_injection = -delta_S_full;

        case 'injection'
            delta_S_injection = delta_S_full;

        otherwise
            error('Unknown game_convention: %s', game_convention);
    end

    Sspec_game = Sspec_stressed + delta_S_injection;

    fprintf('\n============================================================\n');
    fprintf('CASE 3: STRESSED + GAME NONLINEAR AC POWER FLOW\n');
    fprintf('============================================================\n');

    game = solve_acpf_3ph_rect( ...
        Y, slack_idx, Vslack, Sspec_game, stressed.V, solver_opts);

    assert_converged(game, 'stressed + game');
end

%% Summaries

bus_labels = { ...
    '650', '632', '633', '634', '645', '646', '671', ...
    '692', '675', '680', '684', '611', '652' ...
};

phase_labels = {'a', 'b', 'c'};

nominal.summary = summarize_case( ...
    nominal, ns_idx, bus_labels, phase_labels, v_min, v_max);

stressed.summary = summarize_case( ...
    stressed, ns_idx, bus_labels, phase_labels, v_min, v_max);

fprintf('\n============================================================\n');
fprintf('NONLINEAR AC VALIDATION SUMMARY\n');
fprintf('============================================================\n');

print_summary('Nominal', nominal);
print_summary(sprintf('Stress %.2fx', stress_factor), stressed);

if has_game_case
    game.summary = summarize_case( ...
        game, ns_idx, bus_labels, phase_labels, v_min, v_max);
    print_summary('Stress + game', game);
end

%% Plot voltage magnitudes by phase

figure('Name', 'Nonlinear three-phase AC validation');

for phase = 1:3
    subplot(3, 1, phase);

    phase_idx = phase:3:N;

    plot(1:n_bus, nominal.v(phase_idx), '-o', ...
        'DisplayName', 'Nominal');
    hold on;

    plot(1:n_bus, stressed.v(phase_idx), '-s', ...
        'DisplayName', sprintf('Stress %.2fx', stress_factor));

    if has_game_case
        plot(1:n_bus, game.v(phase_idx), '-d', ...
            'DisplayName', 'Stress + game');
    end

    yline(v_min, '--', 'V_{min}');
    yline(v_max, '--', 'V_{max}');

    grid on;
    ylabel(sprintf('|V_%s| [p.u.]', phase_labels{phase}));
    xticks(1:n_bus);
    xticklabels(bus_labels);

    if phase == 1
        title('Nonlinear AC voltage magnitudes');
        legend('Location', 'best');
    end

    if phase == 3
        xlabel('IEEE 13-node feeder bus');
    end
end

%% Save results

results = struct();
results.settings.stress_factor = stress_factor;
results.settings.v_min = v_min;
results.settings.v_max = v_max;
results.settings.game_file = game_file;
results.settings.game_convention = game_convention;

results.nominal = nominal;
results.stressed = stressed;

if has_game_case
    results.game = game;
end

save('nonlinear_ac_validation_results.mat', 'results');

fprintf('\nSaved nonlinear_ac_validation_results.mat\n');

%% Local functions

function assert_converged(result, case_name)
    if ~result.converged
        error([ ...
            'The %s nonlinear AC case did not converge. ', ...
            'Final mismatch = %.3e p.u.' ...
        ], case_name, result.max_mismatch);
    end
end


function summary = summarize_case( ...
    result, ns_idx, bus_labels, phase_labels, v_min, v_max)

    v_ns = result.v(ns_idx);

    [summary.Vmin, local_min_idx] = min(v_ns);
    [summary.Vmax, local_max_idx] = max(v_ns);

    full_min_idx = ns_idx(local_min_idx);
    full_max_idx = ns_idx(local_max_idx);

    [summary.Vmin_bus, summary.Vmin_phase] = ...
        terminal_name(full_min_idx, bus_labels, phase_labels);

    [summary.Vmax_bus, summary.Vmax_phase] = ...
        terminal_name(full_max_idx, bus_labels, phase_labels);

    summary.n_undervoltage = sum(v_ns < v_min);
    summary.n_overvoltage = sum(v_ns > v_max);
    summary.n_voltage_violations = ...
        summary.n_undervoltage + summary.n_overvoltage;

    summary.max_mismatch = result.max_mismatch;
    summary.slack_P_total = sum(real(result.slack_power));
    summary.slack_Q_total = sum(imag(result.slack_power));

    % With a passive series network, sum of calculated injections equals loss.
    summary.P_loss = sum(real(result.S_calc));
    summary.Q_loss = sum(imag(result.S_calc));
end


function [bus_name, phase_name] = terminal_name( ...
    terminal_idx, bus_labels, phase_labels)

    bus_idx = ceil(terminal_idx/3);
    phase_idx = mod(terminal_idx - 1, 3) + 1;

    bus_name = bus_labels{bus_idx};
    phase_name = phase_labels{phase_idx};
end


function print_summary(case_name, result)
    s = result.summary;

    fprintf('\n%s\n', case_name);
    fprintf('  converged             : %d\n', result.converged);
    fprintf('  Newton iterations     : %d\n', result.iterations);
    fprintf('  max mismatch          : %.3e p.u.\n', s.max_mismatch);
    fprintf('  minimum voltage       : %.6f p.u. at bus %s phase %s\n', ...
        s.Vmin, s.Vmin_bus, s.Vmin_phase);
    fprintf('  maximum voltage       : %.6f p.u. at bus %s phase %s\n', ...
        s.Vmax, s.Vmax_bus, s.Vmax_phase);
    fprintf('  voltage violations    : %d\n', s.n_voltage_violations);
    fprintf('  slack P total         : %.6f p.u.\n', s.slack_P_total);
    fprintf('  slack Q total         : %.6f p.u.\n', s.slack_Q_total);
    fprintf('  active-power losses   : %.6f p.u.\n', s.P_loss);
    fprintf('  reactive-power losses : %.6f p.u.\n', s.Q_loss);
end

% in example_ieee13.m or a new script after running both
