%% RUN_IO_MAP_ACCURACY
%
% Quantifies the THIRD row of the accuracy summary that was previously left
% unquantified: how the local input-output map
%
%     v_ns approx v*'' + S_v * Delta_s_ns          (Delta_s_P = 0)
%
% degrades as Delta_s_F moves away from the S'' anchor point, compared
% against a FRESH independent nonlinear AC solve at each perturbed point
% (not the same AC solve used to build the tangent plane).
%
%   e_max(0)     = 0                          by construction (verified, not assumed)
%   e_max(Delta) = max_i |v_pred_i - v_AC_i|   at increasing ||Delta_s_F||_2
%
% Sweep design: for each of several random directions in the 14-dim
% community flexibility space, scale along that direction as a fraction of
% how far the assumed local capability box (ieee13_local_capability.m)
% extends in that direction, so magnitudes range from the anchor (0) out to
% the edge of the assumed capability set (1.0x) -- exactly the range the
% exported A_F_sh/b_F_sh and any downstream dispatch are meant to be used
% over. The previously computed x_feas and x_min (from
% run_feasibility_check.m) are included as named reference points.

clear all
close all
clc

run_ieee13_pipeline;   % fresh S, S', S'', sh, Y, n, der_meta, valid_mask, ...

[A_loc, b_loc, lo, hi] = ieee13_local_capability(der_meta.community); %#ok<ASGLU>
nF = numel(lo);

% Community bus-phase full-vector indices, same order as A_F_sh columns.
idxF_full = [];
for k = 1:numel(der_meta.community)
    c = der_meta.community(k);
    idxF_full = [idxF_full; rw(c.bus, c.phase)]; %#ok<AGROW>
end

idx = sh.idx;
valid_ns = sh.valid_ns;
vstar_ns = Spp.ac.v(idx.nonSlackPhase);   % == tp_ac.v_star(nonSlackPhase), the anchor

%% Build the perturbation sample set

rand('seed', 42); randn('seed', 42);   %#ok<RAND> reproducible sweep

n_dirs = 8;
fractions = [0.05, 0.10, 0.25, 0.50, 0.75, 1.00];

samples = struct('label', {}, 'dF', {});

samples(end+1) = struct('label', 'anchor (Delta=0)', 'dF', zeros(nF,1));

for di = 1:n_dirs
    d = randn(nF,1);
    d = d / norm(d);
    % how far the box extends along +d and -d
    t_pos = min( (hi(d>0)./d(d>0)) );
    if any(d<0); t_pos = min([t_pos; lo(d<0)./d(d<0)]); end
    for fi = 1:numel(fractions)
        dF = fractions(fi) * t_pos * d;
        samples(end+1) = struct('label', sprintf('dir%d @ %.0f%% box', di, 100*fractions(fi)), 'dF', dF); %#ok<AGROW>
    end
end

% Named reference points from the feasibility checkpoint, if available.
if exist(fullfile('data','feasibility_results.mat'), 'file')
    fr = load(fullfile('data','feasibility_results.mat'));
    samples(end+1) = struct('label', 'x_feas (feasibility_lp)', 'dF', fr.feas.x);
    samples(end+1) = struct('label', 'x_min (min_effort_qp)',   'dF', fr.qp.x);
end

fprintf('\n############################################################\n');
fprintf('IO-MAP ACCURACY SWEEP: predicted (tangent at S'''') vs fresh nonlinear AC\n');
fprintf('############################################################\n');
fprintf('%d samples (%d directions x %d fractions + anchor + %d named points)\n', ...
    numel(samples), n_dirs, numel(fractions), numel(samples)-1-n_dirs*numel(fractions));

%% Sweep

results = struct('label',{}, 'norm_dF',{}, 'max_err',{}, 'rms_err',{}, 'converged',{});

for si = 1:numel(samples)
    dF = samples(si).dF;

    p_pert = p_S2; q_pert = q_S2;
    p_pert(idxF_full) = p_pert(idxF_full) + dF(1:nF/2);
    q_pert(idxF_full) = q_pert(idxF_full) + dF(nF/2+1:end);

    % Predicted voltage from the tangent map, no AC solve involved. Build
    % the full non-slack deviation vector directly from p_pert/q_pert
    % relative to the S'' baseline (correct and unambiguous, avoids any
    % column-order bookkeeping error):
    sstar_ns = [Spp.p(idx.nonSlackPhase); Spp.q(idx.nonSlackPhase)];
    s_pert_ns = [p_pert(idx.nonSlackPhase); q_pert(idx.nonSlackPhase)];
    v_pred_ns = sh.S_v*s_pert_ns + sh.m_v;

    % Fresh, independent nonlinear AC solve at the perturbed point.
    Sspec = -(p_pert + 1j*q_pert); Sspec(1:3) = 0;
    Vslack = v_slack .* exp(1j*t_slack);
    ac = solve_acpf_3ph_rect(Y, (1:3)', Vslack, Sspec, Spp.ac.V, ac_opts);

    v_ac_ns = ac.v(idx.nonSlackPhase);
    err = abs(v_pred_ns(valid_ns) - v_ac_ns(valid_ns));

    results(si).label = samples(si).label;
    results(si).norm_dF = norm(s_pert_ns - sstar_ns);
    results(si).max_err = max(err);
    results(si).rms_err = sqrt(mean(err.^2));
    results(si).converged = ac.converged;

    fprintf('  %-24s  ||Delta_s||=%.5f  e_max=%.3e  e_rms=%.3e  (AC converged: %d)\n', ...
        results(si).label, results(si).norm_dF, results(si).max_err, results(si).rms_err, ac.converged);
end

%% Anchor sanity check: error must be exactly (numerically) zero at Delta=0

assert(results(1).max_err < 1e-9, 'run_io_map_accuracy:anchorNotZero', ...
    'Error at the anchor point is %.3e, expected ~0 by construction.', results(1).max_err);
fprintf('\nAnchor check passed: e_max(Delta=0) = %.3e (zero by construction, verified numerically)\n', ...
    results(1).max_err);

%% Plot: error vs perturbation magnitude

norm_dF_all = [results.norm_dF];
max_err_all = [results.max_err];
rms_err_all = [results.rms_err];

figure('visible', 'off');
loglog(norm_dF_all(norm_dF_all>0), max_err_all(norm_dF_all>0), 'ro', 'MarkerFaceColor', 'r'); hold on;
loglog(norm_dF_all(norm_dF_all>0), rms_err_all(norm_dF_all>0), 'bo', 'MarkerFaceColor', 'b');
xlabel('|| Delta s_F ||_2  [pu]');
ylabel('voltage error vs fresh nonlinear AC  [pu]');
title('IO-map (tangent at S'''') accuracy vs distance from anchor');
legend('max error', 'RMS error', 'Location', 'northwest');
grid on;
print(fullfile('data','io_map_accuracy.png'), '-dpng');

save(fullfile('data','io_map_accuracy_results.mat'), 'results', 'samples', '-v7');
fprintf('\nSaved data/io_map_accuracy_results.mat and data/io_map_accuracy.png\n');
