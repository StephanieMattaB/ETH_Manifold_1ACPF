%% RUN_FEASIBILITY_CHECK
%
% Pure feasibility test of the four-player energy-community flexibility
% problem at the corrected S'' interface -- NO game objectives, NO Nash
% equilibrium, NO Q/q cost structure, NO Nabetani/PWA machinery. Only:
%
%   F = { x_F : A_F_sh*x_F <= b_F_sh,  A_loc*x_F <= b_loc }
%
% Checkpoint order (per the agreed structure):
%   1. re-derive the corrected S'' interface (run_ieee13_pipeline.m)
%   2. confirm x_F = 0 is infeasible w.r.t. the AC-derived violated rows
%   3. solve the feasibility LP: is F empty or not?
%   4. if non-empty, solve the minimum-effort QP: min (1/2)||x_F||^2 s.t. F
%   5. validate the minimum-effort dispatch in the NONLINEAR AC model
%
% IMPORTANT SCOPING NOTE: this script is a DOWNSTREAM DIAGNOSTIC of the DSO
% output (A_F_sh, b_F_sh from derive_shared_constraint.m), not a
% precondition for it. A_F_sh/b_F_sh are already fully valid as the LOCAL
% voltage-security constraint of the community's controllable variables
% around the selected, AC-validated S'' operating point -- that is exactly
% the theoretical meaning of a tangent-plane constraint, and nothing here
% changes it. What steps 4-5 test is a SEPARATE question: whether one
% particular (minimum-effort) corrective dispatch, evaluated far enough
% from the S'' linearization point to matter, still holds up when re-solved
% in the true nonlinear AC model. A negative answer to that question is a
% statement about the RANGE OF VALIDITY of the local constraint under a
% large correction -- it does not retract S'' as the adopted operating
% point, nor invalidate A_F_sh/b_F_sh as the DSO-to-community interface.

clear all
close all
clc

run_ieee13_pipeline;   % re-derives S, S', S'', sh (corrected interface), Y, n,
                        % der_meta, valid_mask, vmin, vmax, ac_opts fresh

fprintf('\n\n############################################################\n');
fprintf('FEASIBILITY CHECK: F = {x_F : A_F_sh x_F <= b_F_sh, A_loc x_F <= b_loc}\n');
fprintf('############################################################\n');

%% Step 2: x_F = 0 infeasibility, confirmed directly from the (corrected) RHS

n_neg_upper = nnz(sh.b_F_sh(1:sh.nPhys) < 0);
n_neg_lower = nnz(sh.b_F_sh(sh.nPhys+1:end) < 0);
x0_infeasible = (n_neg_upper + n_neg_lower) > 0;

fprintf('\n-- Step 2: x_F = 0 feasibility (shared constraint alone) --\n');
fprintf('  b_F_sh negative rows: %d upper, %d lower (of %d each)\n', ...
    n_neg_upper, n_neg_lower, sh.nPhys);
fprintf('  x_F = 0 is %s w.r.t. the shared voltage-security constraint\n', ...
    ternary(x0_infeasible, 'INFEASIBLE', 'feasible'));
fprintf('  (must match AC ground truth: %d over / %d under non-slack violations found earlier)\n', ...
    Spp.rep_ac_ns.n_overvoltage, Spp.rep_ac_ns.n_undervoltage);
assert(n_neg_upper == Spp.rep_ac_ns.n_overvoltage && n_neg_lower == Spp.rep_ac_ns.n_undervoltage, ...
    'run_feasibility_check:mismatch', 'x_F=0 infeasibility pattern does not match AC ground truth.');

%% Step 3: local capability constraints + feasibility LP

[A_loc, b_loc, lo, hi] = ieee13_local_capability(der_meta.community);

fprintf('\n-- Step 3: local capability constraints --\n');
fprintf('  A_loc: %d x %d (box on all 14 [p,q] community variables)\n', size(A_loc,1), size(A_loc,2));
fprintf('  bounds (pu): p in [-0.20,+1.00]xP_DER_i, q in [-0.30,+0.30]xP_DER_i (see ieee13_local_capability.m)\n');

feas = feasibility_lp(sh.A_F_sh, sh.b_F_sh, A_loc, b_loc);

fprintf('\n-- Feasibility LP result --\n');
fprintf('  solver status          : %d (5 = optimal)\n', feas.status);
fprintf('  min worst-row slack t* : %.6f pu\n', feas.t);
fprintf('  F is %s\n', ternary(feas.feasible, 'NON-EMPTY (feasible)', 'EMPTY (infeasible)'));

if ~feas.feasible
    fprintf('\n>>> F is EMPTY: the assumed local capability set cannot restore all shared\n');
    fprintf('>>> voltage-security rows simultaneously. This means S'''' is too stressed for\n');
    fprintf('>>> these four participants under the assumed capability -- reconsider lambda,\n');
    fprintf('>>> DER sizing, or the capability assumption before building the game.\n');
    return;
end

fprintf('  worst shared-row slack at x_feas: %.6f pu\n', min(feas.slack_sh));
fprintf('  worst local-row slack at x_feas : %.6f pu\n', min(feas.slack_loc));
fprintf('  ||x_feas||_2                     : %.6f pu\n', norm(feas.x));

n_resolved_sh = nnz(sh.b_F_sh < 0 & feas.slack_sh >= -1e-9);
fprintf('  originally-violated shared rows now satisfied: %d / %d\n', n_resolved_sh, n_neg_upper+n_neg_lower);

%% Step 4: minimum-effort corrective dispatch

qp = min_effort_qp(sh.A_F_sh, sh.b_F_sh, A_loc, b_loc);

fprintf('\n-- Step 4: minimum-effort dispatch  min (1/2)||x_F||^2  s.t. F --\n');
fprintf('  ADMM iterations         : %d (converged: %d)\n', qp.iterations, qp.converged);
fprintf('  max constraint violation: %.3e (must be <= 0)\n', qp.max_constraint_violation);
fprintf('  ||x_min||_2              : %.6f pu   (vs ||x_feas||_2 = %.6f pu)\n', qp.norm2, norm(feas.x));

assert(qp.feasible, 'run_feasibility_check:qpInfeasible', ...
    'Minimum-effort QP solution violates constraints -- refusing to proceed to AC validation.');

slack_sh_min = sh.b_F_sh - sh.A_F_sh*qp.x;
n_resolved_sh_min = nnz(sh.b_F_sh < 0 & slack_sh_min >= -1e-9);
fprintf('  originally-violated shared rows now satisfied (min-effort): %d / %d\n', ...
    n_resolved_sh_min, n_neg_upper+n_neg_lower);

fprintf('\n  Per-player minimum-effort dispatch (pu, community column order):\n');
layout = sh.layout;   % single source of truth (community_column_layout.m), NOT re-derived here
for i = 1:layout.n
    fprintf('    %-6s : %+.5f  (bounds [%+.4f, %+.4f])\n', layout.labels{i}, qp.x(i), lo(i), hi(i));
end

%% Step 5: AC-validate the minimum-effort dispatch

fprintf('\n-- Step 5: nonlinear-AC validation of x_min --\n');

% Map x_min (14, layout column order) back onto the full 3n-vector net-
% consumption specification, s'' + x_min: loop column-by-column using
% layout.full_idx/layout.is_q rather than assuming any block structure --
% this is exactly the mapping that was previously wrong (a flat
% p-block-then-q-block split that did not match A_F_sh's actual per-player
% column order), corrupting 9 of the 14 columns' physical meaning.
p_corrected = p_S2;
q_corrected = q_S2;
for i = 1:layout.n
    if layout.is_q(i)
        q_corrected(layout.full_idx(i)) = q_corrected(layout.full_idx(i)) + qp.x(i);
    else
        p_corrected(layout.full_idx(i)) = p_corrected(layout.full_idx(i)) + qp.x(i);
    end
end

corrected = run_stage('S'''' + minimum-effort correction', Y, n, v_slack, t_slack, ...
    p_corrected, q_corrected, valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, Spp.ac.V);

fprintf('\n-- Before/after (nonlinear AC, non-slack physical terminals) --\n');
fprintf('  S'''' baseline   : Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    Spp.rep_ac_ns.Vmin, Spp.rep_ac_ns.Vmax, Spp.rep_ac_ns.n_violations, Spp.rep_ac_ns.n_physical);
fprintf('  after correction: Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    corrected.rep_ac_ns.Vmin, corrected.rep_ac_ns.Vmax, corrected.rep_ac_ns.n_violations, corrected.rep_ac_ns.n_physical);

%% Consistency printout (per the requested cross-check): same x_min, both
%  the LINEAR prediction and a FRESH independent AC solve, with the row
%  attaining each minimum -- must satisfy V_min_lin >= vmin (by feasibility
%  of x_min) and |V_min_lin - V_min_AC_fresh| <= e_max (by definition of
%  e_max), so V_min_AC_fresh >= vmin - e_max necessarily.

s_corrected_ns = [p_corrected(sh.idx.nonSlackPhase); q_corrected(sh.idx.nonSlackPhase)];
v_pred_ns = sh.S_v*s_corrected_ns + sh.m_v;
v_pred_phys = v_pred_ns(sh.valid_ns);
[v_pred_min, i_pred_min] = min(v_pred_phys);

v_ac_fresh_ns = corrected.ac.v(sh.idx.nonSlackPhase);
v_ac_fresh_phys = v_ac_fresh_ns(sh.valid_ns);
[v_ac_min, i_ac_min] = min(v_ac_fresh_phys);

err_phys = abs(v_pred_phys - v_ac_fresh_phys);
[e_max, i_err_max] = max(err_phys);

phys_full_idx = sh.idx.nonSlackPhase(sh.valid_ns);
term = @(k) sprintf('%s.%s', bus_labels{ceil(phys_full_idx(k)/3)}, phase_labels{mod(phys_full_idx(k)-1,3)+1});

fprintf('\n-- Consistency check (same x_min, linear prediction vs fresh AC) --\n');
fprintf('  ||x_min||_2                    = %.6f pu\n', norm(qp.x));
fprintf('  min(v_lin(x_min))    = %.6f pu at %s   (must be >= vmin=%.2f since x_min in F)\n', ...
    v_pred_min, term(i_pred_min), vmin);
fprintf('  min(v_AC_fresh(x_min)) = %.6f pu at %s\n', v_ac_min, term(i_ac_min));
fprintf('  max|v_lin - v_AC_fresh|         = %.3e pu at %s\n', e_max, term(i_err_max));
fprintf('  implied lower bound on v_AC_fresh: vmin - e_max = %.6f pu\n', vmin - e_max);
assert(v_pred_min >= vmin - 1e-9, 'run_feasibility_check:linearInfeasible', ...
    'x_min''s own linear prediction violates vmin -- x_min is not actually in F.');
assert(v_ac_min >= vmin - e_max - 1e-9, 'run_feasibility_check:inconsistentWithErrorBound', ...
    'v_AC_fresh(x_min) = %.6f is below vmin-e_max = %.6f -- contradicts the error bound, mapping bug likely.', ...
    v_ac_min, vmin - e_max);
fprintf('  Consistency check PASSED: no contradiction between linear feasibility and the AC result.\n');

if corrected.rep_ac_ns.n_violations == 0
    fprintf('\n>>> CERTIFIED: the minimum-effort linear correction, applied to the true nonlinear\n');
    fprintf('>>> AC model, restores 0.95 <= |V| <= 1.05 at every non-slack physical terminal.\n');
else
    fprintf('\n>>> NOT CERTIFIED: %d violations remain in the nonlinear AC model after applying\n', ...
        corrected.rep_ac_ns.n_violations);
    fprintf('>>> the linear-tangent-optimal correction -- linear feasibility does NOT guarantee\n');
    fprintf('>>> nonlinear recoverability at this stress level. Remaining violated terminals:\n');
    fprintf('    %s\n', strjoin(corrected.rep_ac_ns.violated_terminals, ', '));
end

save(fullfile('data', 'feasibility_results.mat'), 'feas', 'qp', 'A_loc', 'b_loc', 'lo', 'hi', ...
    'corrected', 'p_corrected', 'q_corrected', '-v7');
fprintf('\nSaved feasibility/min-effort/AC-validation results to data/feasibility_results.mat\n');
