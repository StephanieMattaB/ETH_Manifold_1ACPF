%% VALIDATE_FAIRNESS_GNE
%
% Nonlinear AC validation of the externally-supplied minimax-fair GNE
% x*_fair (from the Julia constrained-game solve, region CR_1,
% fairness_selected_x.csv), following exactly the same procedure as
% Step 5 of run_feasibility_check.m: map x* onto the S'' network via
% community_column_layout.m, re-solve the TRUE nonlinear three-phase AC
% power flow at the corrected dispatch (fresh solve, not reusing the
% tangent-plane prediction), and report before/after voltage security
% plus the linear-vs-AC consistency check.
%
% This script does NOT recompute x*; it takes it as given, verified
% input from the Julia side and asks only: does it hold up in the
% original nonlinear network?

clear all
close all
clc

run_ieee13_pipeline;   % re-derives S, S', S'', sh (corrected interface), Y, n,
                        % der_meta, valid_mask, vmin, vmax, ac_opts, layout,
                        % p_S2/q_S2 fresh

layout = sh.layout;    % single source of truth (community_column_layout.m)

fprintf('\n\n############################################################\n');
fprintf('NONLINEAR AC VALIDATION OF EXTERNALLY-SUPPLIED FAIRNESS GNE x*_fair\n');
fprintf('Source: fairness_selected_x.csv, region CR_1\n');
fprintf('############################################################\n');

%% x*_fair, in EXACT community_column_layout.m order
% [p645b,p645c,q645b,q645c, p611c,q611c, p652a,q652a, p671a,p671b,p671c,q671a,q671b,q671c]

x_star = [ ...
    -0.001127; ...  % p645b
     0.002666; ...  % p645c
     0.010958; ...  % q645b
    -0.018257; ...  % q645c
    -0.008425; ...  % p611c
    -0.021260; ...  % q611c
    -0.013423; ...  % p652a
    -0.018194; ...  % q652a
    -0.017823; ...  % p671a
     0.001390; ...  % p671b
    -0.001349; ...  % p671c
    -0.018800; ...  % q671a
     0.010285; ...  % q671b
    -0.026691  ...  % q671c
];

assert(numel(x_star) == layout.n, 'validate_fairness_gne:sizeMismatch', ...
    'x_star has %d entries, expected %d (layout.n).', numel(x_star), layout.n);

fprintf('\n-- Input check --\n');
fprintf('  ||x_star||_2 (as supplied)        : %.6f pu\n', norm(x_star));
fprintf('  ||x_star||_2 (recomputed here)    : %.6f pu  (must match)\n', norm(x_star));
for i = 1:layout.n
    fprintf('    %-6s = %+.6f pu  (%+.3f kW/kVAr)\n', layout.labels{i}, x_star(i), x_star(i)*5000);
end

%% Step 1: local capability box check

[A_loc, b_loc, lo, hi] = ieee13_local_capability(der_meta.community);
loc_violation = max([lo - x_star; x_star - hi]);
fprintf('\n-- Step 1: local capability box (ieee13_local_capability.m) --\n');
fprintf('  max violation of [lo,hi] box       : %.3e pu  (<=0 means inside the box)\n', loc_violation);
for i = 1:layout.n
    if x_star(i) < lo(i) - 1e-9 || x_star(i) > hi(i) + 1e-9
        fprintf('    OUT OF BOX: %-6s = %+.6f  not in [%+.4f, %+.4f]\n', ...
            layout.labels{i}, x_star(i), lo(i), hi(i));
    end
end
if loc_violation <= 1e-9
    fprintf('  x_star is INSIDE the local capability box.\n');
else
    fprintf('  >>> x_star VIOLATES the local capability box by %.3e pu.\n', loc_violation);
end

%% Step 2: shared (linear) constraint check, at the corrected S'' interface

shared_slack = sh.b_F_sh - sh.A_F_sh * x_star;
fprintf('\n-- Step 2: shared voltage-security constraint (A_F_sh x <= b_F_sh) --\n');
fprintf('  worst-row slack (min over 58 rows) : %.6f pu  (>=0 means linearly feasible)\n', min(shared_slack));
if min(shared_slack) >= -1e-9
    fprintf('  x_star is LINEARLY FEASIBLE w.r.t. the shared constraint.\n');
else
    fprintf('  >>> x_star VIOLATES the shared constraint by %.3e pu (linear prediction).\n', -min(shared_slack));
end

%% Step 3: map x_star onto the full network (S'' + x_star), exactly as
%  run_feasibility_check.m Step 5 does -- loop by layout.is_q/full_idx,
%  never a flat p-then-q split.

p_corrected = p_S2;
q_corrected = q_S2;
for i = 1:layout.n
    if layout.is_q(i)
        q_corrected(layout.full_idx(i)) = q_corrected(layout.full_idx(i)) + x_star(i);
    else
        p_corrected(layout.full_idx(i)) = p_corrected(layout.full_idx(i)) + x_star(i);
    end
end

%% Step 4: fresh, independent nonlinear AC solve at the corrected dispatch

corrected = run_stage('S'''' + fairness-GNE correction (x*_fair)', Y, n, v_slack, t_slack, ...
    p_corrected, q_corrected, valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, Spp.ac.V);

fprintf('\n-- Before/after (nonlinear AC, non-slack physical terminals) --\n');
fprintf('  S'''' baseline    : Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    Spp.rep_ac_ns.Vmin, Spp.rep_ac_ns.Vmax, Spp.rep_ac_ns.n_violations, Spp.rep_ac_ns.n_physical);
fprintf('  after x*_fair    : Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    corrected.rep_ac_ns.Vmin, corrected.rep_ac_ns.Vmax, corrected.rep_ac_ns.n_violations, corrected.rep_ac_ns.n_physical);

%% Step 5: consistency check -- same x_star, linear prediction vs fresh AC

s_corrected_ns = [p_corrected(sh.idx.nonSlackPhase); q_corrected(sh.idx.nonSlackPhase)];
v_pred_ns = sh.S_v*s_corrected_ns + sh.m_v;
v_pred_phys = v_pred_ns(sh.valid_ns);
[v_pred_min, i_pred_min] = min(v_pred_phys);
[v_pred_max, i_pred_max] = max(v_pred_phys);

v_ac_fresh_ns = corrected.ac.v(sh.idx.nonSlackPhase);
v_ac_fresh_phys = v_ac_fresh_ns(sh.valid_ns);
[v_ac_min, i_ac_min] = min(v_ac_fresh_phys);
[v_ac_max, i_ac_max] = max(v_ac_fresh_phys);

err_phys = abs(v_pred_phys - v_ac_fresh_phys);
[e_max, i_err_max] = max(err_phys);
e_rms = sqrt(mean(err_phys.^2));

phys_full_idx = sh.idx.nonSlackPhase(sh.valid_ns);
term = @(k) sprintf('%s.%s', bus_labels{ceil(phys_full_idx(k)/3)}, phase_labels{mod(phys_full_idx(k)-1,3)+1});

fprintf('\n-- Consistency check (same x*_fair, linear prediction vs fresh AC) --\n');
fprintf('  min(v_lin(x*))     = %.6f pu at %s\n', v_pred_min, term(i_pred_min));
fprintf('  min(v_AC_fresh(x*))= %.6f pu at %s\n', v_ac_min, term(i_ac_min));
fprintf('  max(v_lin(x*))     = %.6f pu at %s\n', v_pred_max, term(i_pred_max));
fprintf('  max(v_AC_fresh(x*))= %.6f pu at %s\n', v_ac_max, term(i_ac_max));
fprintf('  max|v_lin - v_AC_fresh| = %.3e pu at %s\n', e_max, term(i_err_max));
fprintf('  RMS|v_lin - v_AC_fresh| = %.3e pu\n', e_rms);
fprintf('  implied bound on v_AC_fresh: [vmin-e_max, vmax+e_max] = [%.6f, %.6f]\n', ...
    vmin - e_max, vmax + e_max);

%% Step 6: final certification verdict

fprintf('\n-- Final verdict --\n');
if corrected.rep_ac_ns.n_violations == 0
    fprintf('>>> CERTIFIED: x*_fair, applied to the TRUE nonlinear AC model, restores\n');
    fprintf('>>> 0.95 <= |V| <= 1.05 at every one of the %d non-slack physical terminals.\n', ...
        corrected.rep_ac_ns.n_physical);
else
    fprintf('>>> NOT FULLY CERTIFIED: %d/%d violations remain in the nonlinear AC model\n', ...
        corrected.rep_ac_ns.n_violations, corrected.rep_ac_ns.n_physical);
    fprintf('>>> after applying x*_fair. Remaining violated terminals:\n');
    fprintf('    %s\n', strjoin(corrected.rep_ac_ns.violated_terminals, ', '));
end

save(fullfile('data', 'fairness_gne_validation.mat'), 'x_star', 'corrected', ...
    'p_corrected', 'q_corrected', 'v_pred_phys', 'v_ac_fresh_phys', 'e_max', 'e_rms', '-v7');
fprintf('\nSaved data/fairness_gne_validation.mat\n');
