%% IEEE13_AC_VALIDATION_FINAL
%
% FINAL VALIDATION (replay-only). Re-uses the existing nonlinear unbalanced
% AC pipeline (run_ieee13_pipeline.m / run_stage.m / solve_acpf_3ph_rect.m)
% to independently re-solve three FINAL Julia game equilibria at the
% stressed S'' operating point, and quantifies the error of the linearized
% DSO tangent model (sh.S_v / sh.A_F_sh / sh.vstar_phys / sh.b_F_sh, built
% by derive_shared_constraint.m at S'''s AC-converged point) against the
% true nonlinear AC solution.
%
% This script does NOT re-derive the feeder, DOES NOT build a new
% optimization problem, and DOES NOT change kappa/capability bounds/voltage
% limits/load stress/DER locations/capacitor states/sign conventions. It
% only:
%   1. rebuilds S'' exactly as run_ieee13_pipeline.m does (unchanged call);
%   2. adds each equilibrium's 14 Delta-p/Delta-q deviations onto the S''
%      community bus-phases, using EXACTLY the column mapping already
%      established as the single source of truth (community_column_layout.m,
%      reused verbatim from run_feasibility_check.m's pattern);
%   3. re-solves the FULL nonlinear AC model (run_stage.m, unchanged);
%   4. evaluates the SAME tangent linear model already exported to Julia
%      (sh.A_F_sh/sh.b_F_sh/sh.vstar_phys, i.e. v_lin(x) = vstar_phys +
%      A_F_sh(1:nPhys,:)*x) at the same 29 physical non-slack bus-phase
%      terminals, in the identical ordering;
%   5. reports the linearization error and writes the CSV/summary outputs.
%
% CRITICAL SIGN CONVENTION (verified against the existing code, not assumed):
%   ieee13.m's p,q are NET CONSUMPTION (= load - generation), positive =
%   drawn from the network (see solve_acpf_3ph_rect.m: Sspec = -(p+jq)).
%   ieee13_der.m / ieee13_congest.m subtract DER capacity from net
%   consumption: P_net = P_load - P_DER, Q_net = Q_load - Q_DER (Q_DER,0=0,
%   unity-PF baseline). The game's convention states positive Delta p/Delta q
%   = increased net consumption; since P_net,new = P_net,old + Delta p and
%   P_net,old = P_load - P_DER,0, holding P_load fixed gives
%       P_DER,new = P_DER,0 - Delta p,   Q_DER,new = Q_DER,0 - Delta q = -Delta q.
%   Equivalently and MINIMALLY: simply add Delta p/Delta q directly onto the
%   S'' net-consumption vectors p_S2/q_S2 at the corresponding bus-phase
%   (this is exactly what run_feasibility_check.m already does for x_min).
%   Both are algebraically identical; the P_DER,new/Q_DER,new figures below
%   are reported purely as a sign/mapping sanity check.

clear all
close all
clc

%% ------------------------------------------------------------------------
%  Step 0: rebuild S'' exactly as the existing pipeline does (unchanged).
%  NOTE: run_ieee13_pipeline.m itself starts with `clear all`, so nothing
%  can be set in this workspace before this call -- matches the pattern
%  already used by run_feasibility_check.m.
%  ------------------------------------------------------------------------

run_ieee13_pipeline;   % -> Y, n, v_slack, t_slack, valid_mask, bus_labels,
                        %    phase_labels, vmin, vmax, ac_opts, der_meta,
                        %    p_S2, q_S2, Spp (S'' AC stage), sh (tangent /
                        %    shared-constraint interface at S'''s AC point),
                        %    lambda

results_dir = fullfile('results', 'ac_validation_final');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

layout = sh.layout;    % single source of truth (community_column_layout.m)

fprintf('\n\n############################################################\n');
fprintf('FINAL AC REPLAY: S'''' baseline + 3 final Julia equilibria\n');
fprintf('############################################################\n');

%% ------------------------------------------------------------------------
%  Step 0b: baseline reproduction check (task: do not proceed past this
%  without confirming S'''' reproduces the expected stressed nonlinear point)
%  ------------------------------------------------------------------------

expected.Vmin = 0.905074; expected.Vmax = 1.06067;
expected.n_under = 14; expected.n_over = 4; expected.n_total = 18;

got.Vmin = Spp.rep_ac_ns.Vmin; got.Vmax = Spp.rep_ac_ns.Vmax;
got.n_under = Spp.rep_ac_ns.n_undervoltage; got.n_over = Spp.rep_ac_ns.n_overvoltage;
got.n_total = Spp.rep_ac_ns.n_violations;

tol_v = 1e-3;
baseline_ok = abs(got.Vmin - expected.Vmin) < tol_v && ...
              abs(got.Vmax - expected.Vmax) < tol_v && ...
              got.n_under == expected.n_under && got.n_over == expected.n_over && ...
              got.n_total == expected.n_total;

fprintf('\n-- S'''' baseline reproduction check --\n');
fprintf('  expected: Vmin=%.6f Vmax=%.6f under=%d over=%d total=%d/29\n', ...
    expected.Vmin, expected.Vmax, expected.n_under, expected.n_over, expected.n_total);
fprintf('  got     : Vmin=%.6f Vmax=%.6f under=%d over=%d total=%d/29\n', ...
    got.Vmin, got.Vmax, got.n_under, got.n_over, got.n_total);
fprintf('  baseline reproduced = %s\n', ternary(baseline_ok, 'YES', 'NO'));

if ~baseline_ok
    error('ieee13_ac_validation_final:baselineMismatch', ...
        ['S'''' nonlinear AC baseline does NOT reproduce the expected stressed operating point.\n' ...
         'STOPPING per instructions -- diagnose before proceeding (do not hard-code the expected values).']);
end

%% ------------------------------------------------------------------------
%  Step 0c: physical bus-phase terminal ordering (29 non-slack terminals),
%  identical to the ordering already used by sh.A_F_sh / sh.vstar_phys.
%  ------------------------------------------------------------------------

phys_full_idx = sh.idx.nonSlackPhase(sh.valid_ns);   % 29x1, full 3n indexing
nPhys = sh.nPhys;
term = @(k) sprintf('%s.%s', bus_labels{ceil(phys_full_idx(k)/3)}, phase_labels{mod(phys_full_idx(k)-1,3)+1});
term_bus = @(k) bus_labels{ceil(phys_full_idx(k)/3)};
term_phase = @(k) phase_labels{mod(phys_full_idx(k)-1,3)+1};

fprintf('\n-- Active linear shared constraint rows (from Julia POINT A) --\n');
for r = [49, 51]
    if r <= nPhys
        k = r; kind = 'upper (vmax)';
    else
        k = r - nPhys; kind = 'lower (vmin)';
    end
    fprintf('  row %d = %s constraint at physical terminal #%d = %s\n', r, kind, k, term(k));
end

%% ------------------------------------------------------------------------
%  Step 1: define the FINAL equilibrium points (verbatim from the Julia
%  results; column order keyed by label, NOT by position, so it is immune
%  to any accidental reordering).
%  ------------------------------------------------------------------------

points = struct('label', {}, 'shortname', {}, 'x_map', {}, 'julia', {});

% ---- POINT A: constrained vGNE ------------------------------------------
mA = containers.Map();
mA('p645b') =  9.797244e-04;  mA('p645c') = -6.438873e-04;
mA('q645b') =  5.073877e-03;  mA('q645c') = -1.242916e-02;
mA('p611c') = -1.622434e-03;  mA('q611c') = -2.594287e-02;
mA('p652a') = -1.685022e-02;  mA('q652a') = -1.789296e-02;
mA('p671a') = -1.689925e-02;  mA('p671b') =  1.564403e-03;  mA('p671c') = -1.634707e-03;
mA('q671a') = -1.743309e-02;  mA('q671b') =  9.698941e-03;  mA('q671c') = -2.590583e-02;

points(end+1) = struct('label', 'constrained vGNE', 'shortname', 'vGNE', ...
    'x_map', mA, 'julia', struct('E2', 1.410577e-03, 'Vmin_lin', 0.95000, ...
    'Vmax_lin', 1.02208, 'viol_lin', 0, 'binding_rows', [49 51], ...
    'max_local_util', 0.929335));

% ---- POINT B: selected compromise (alpha = 0.075) ------------------------
% Source: results/ieee13_constrained_game/selection_20260812_205251/
%         selected_equilibria_full.jls, field obj.compromise.x (Julia
%         Serialization.deserialize), read and reported verbatim -- NOT
%         reconstructed from alpha/printed costs. Provenance cross-checked
%         against alpha=0.075, region_id=1, excess_effort=0.1298%,
%         reference_burden_relief=7.3207%, max_local_utilization=0.960466,
%         shared_feasible=true, local_feasible=true (all match the task spec).
mB = containers.Map();
mB('p645b') =  1.0037729157e-03;  mB('p645c') = -6.5312771706e-04;
mB('q645b') =  5.2556971262e-03;  mB('q645c') = -1.2845270832e-02;
mB('p611c') = -1.6500943409e-03;  mB('q611c') = -2.6811910877e-02;
mB('p652a') = -1.7447459733e-02;  mB('q652a') = -1.8560431509e-02;
mB('p671a') = -1.6277594987e-02;  mB('p671b') =  1.4879060343e-03;  mB('p671c') = -1.5472663273e-03;
mB('q671a') = -1.6822002252e-02;  mB('q671b') =  9.3462101463e-03;  mB('q671c') = -2.4906893342e-02;

points(end+1) = struct('label', 'selected compromise', 'shortname', 'compromise', ...
    'x_map', mB, 'julia', struct('E2', NaN, 'Vmin_lin', NaN, 'Vmax_lin', NaN, ...
    'viol_lin', 0, 'binding_rows', [], 'max_local_util', 0.960466));

% ---- POINT C: high-relief selected point (alpha ~= 1) ---------------------
% Source: same file, field obj.burden_endpoint.x. Provenance cross-checked
% against alpha=1.0, region_id=1, excess_effort=0.6395%,
% reference_burden_relief=15.5987%, max_local_utilization=0.988254,
% shared_feasible=true, local_feasible=true (all match the task spec).
mC = containers.Map();
mC('p645b') =  9.5246670559e-04;  mC('p645c') = -5.5911690717e-04;
mC('q645b') =  5.5160292277e-03;  mC('q645c') = -1.3214051860e-02;
mC('p611c') = -1.4532564596e-03;  mC('q611c') = -2.7587610719e-02;
mC('p652a') = -1.8251830033e-02;  mC('q652a') = -1.9721334832e-02;
mC('p671a') = -1.5464832034e-02;  mC('p671b') =  1.5610361084e-03;  mC('p671c') = -1.6824530526e-03;
mC('q671a') = -1.5747059657e-02;  mC('q671b') =  8.8479793978e-03;  mC('q671c') = -2.4020896386e-02;

points(end+1) = struct('label', 'high-relief selected point', 'shortname', 'high_relief', ...
    'x_map', mC, 'julia', struct('E2', NaN, 'Vmin_lin', NaN, 'Vmax_lin', NaN, ...
    'viol_lin', 0, 'binding_rows', [], 'max_local_util', 0.988254));

if numel(points) < 3
    warning('ieee13_ac_validation_final:missingPoints', ...
        'Only %d/3 final equilibrium points are populated.', numel(points));
end

%% ------------------------------------------------------------------------
%  Step 2: baseline (S'') row for the summary/CSV outputs
%  ------------------------------------------------------------------------

v_ac_S2   = Spp.ac.v(phys_full_idx);
v_lin_S2  = sh.vstar_phys;    % tangent is anchored AT S'' AC point: identical by construction

summary_rows = {};
busphase_rows = {};
injection_rows = {};

summary_rows(end+1,:) = make_summary_row('S'''' stressed (tangent anchor == AC, by construction)', ...
    v_lin_S2, v_ac_S2, vmin, vmax);
busphase_rows = [busphase_rows; make_busphase_rows('S'''' stressed', v_lin_S2, v_ac_S2, ...
    phys_full_idx, bus_labels, phase_labels, vmin, vmax)];

%% ------------------------------------------------------------------------
%  Step 3: for each populated point -- injection mapping, AC replay, linear
%  prediction, error quantification. p_S2/q_S2 (the S'' baseline) are never
%  mutated in place, so deviations are NOT cumulative across points.
%  ------------------------------------------------------------------------

replayed = struct('label', {}, 'x', {}, 'ac', {}, 'v_lin', {}, 'v_ac', {}, 'err', {}, ...
    'ac_feasible', {}, 'case_class', {}, 'worst_viol_mag', {}, 'margin_reco', {});

for ip = 1:numel(points)
    pt = points(ip);
    fprintf('\n============================================================\n');
    fprintf('POINT: %s\n', pt.label);
    fprintf('============================================================\n');

    x = zeros(layout.n, 1);
    for i = 1:layout.n
        x(i) = pt.x_map(layout.labels{i});
    end

    % ---- Task 3: injection mapping table (sign/mapping sanity check) -----
    fprintf('\nlabel        Delta pu        Delta physical     P/Q_DER,new\n');
    fprintf('-----------------------------------------------------------------------\n');
    for i = 1:layout.n
        lab = layout.labels{i};
        k = layout.player_of_col(i);
        c = der_meta.community(k);
        pidx = find(c.phase == layout.phase(i), 1);
        P_DER0_kW = c.p_kW(pidx);
        Q_DER0_kVAr = 0;   % unity-PF baseline, verified in ieee13_der.m
        Sbase_kVA = 5e3;   % 5e6 VA / 1e3

        if layout.is_q(i)
            dQ_kVAr = x(i) * Sbase_kVA;
            Q_DER_new = Q_DER0_kVAr - dQ_kVAr;
            fprintf('%-10s   %+.6e   %+9.4f kVAr   Q_DER,new = %+9.4f kVAr\n', ...
                lab, x(i), dQ_kVAr, Q_DER_new);
            injection_rows(end+1,:) = {pt.shortname, c.bus, phase_labels{layout.phase(i)}, ...
                'Q', x(i), dQ_kVAr, Q_DER0_kVAr, Q_DER_new}; %#ok<AGROW>
        else
            dP_kW = x(i) * Sbase_kVA;
            P_DER_new = P_DER0_kW - dP_kW;
            fprintf('%-10s   %+.6e   %+9.4f kW     P_DER,new = %+9.4f kW\n', ...
                lab, x(i), dP_kW, P_DER_new);
            injection_rows(end+1,:) = {pt.shortname, c.bus, phase_labels{layout.phase(i)}, ...
                'P', x(i), dP_kW, P_DER0_kW, P_DER_new}; %#ok<AGROW>
        end
    end

    % ---- Task 2: apply deviations to a FRESH copy of the S'' baseline ----
    p_pt = p_S2;   % restore-from-original: always start from the untouched
    q_pt = q_S2;   % S'' vectors, never from a previous point's result
    for i = 1:layout.n
        if layout.is_q(i)
            q_pt(layout.full_idx(i)) = q_pt(layout.full_idx(i)) + x(i);
        else
            p_pt(layout.full_idx(i)) = p_pt(layout.full_idx(i)) + x(i);
        end
    end

    % ---- Task 2: full nonlinear AC re-solve (existing run_stage.m) -------
    stage = run_stage(sprintf('S'''' + %s', pt.label), Y, n, v_slack, t_slack, ...
        p_pt, q_pt, valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, Spp.ac.V);

    % ---- Task 4: linear (DSO tangent) prediction, SAME 29 terminals ------
    v_lin_phys = sh.vstar_phys + sh.A_F_sh(1:nPhys, :) * x;
    v_ac_phys  = stage.ac.v(phys_full_idx);
    err = v_lin_phys - v_ac_phys;

    fprintf('\n-- Linear (DSO tangent) vs nonlinear AC, %d physical terminals --\n', nPhys);
    fprintf('  max|error| = %.6e pu   RMS = %.6e pu   mean|error| = %.6e pu\n', ...
        max(abs(err)), sqrt(mean(err.^2)), mean(abs(err)));
    [~, iw] = max(abs(err));
    fprintf('  worst terminal: %s  (v_lin=%.6f, v_AC=%.6f)\n', term(iw), v_lin_phys(iw), v_ac_phys(iw));

    fprintf('\n  cross-check vs Julia-reported linear figures:\n');
    fprintf('    MATLAB v_lin  Vmin=%.6f Vmax=%.6f   |  Julia linear Vmin=%.5f Vmax=%.5f\n', ...
        min(v_lin_phys), max(v_lin_phys), pt.julia.Vmin_lin, pt.julia.Vmax_lin);

    tol = 1e-6;
    n_under_lin = nnz(v_lin_phys < vmin); n_over_lin = nnz(v_lin_phys > vmax);
    n_under_ac  = nnz(v_ac_phys  < vmin); n_over_ac  = nnz(v_ac_phys  > vmax);
    n_under_lin_tol = nnz(v_lin_phys < vmin - tol); n_over_lin_tol = nnz(v_lin_phys > vmax + tol);
    n_under_ac_tol  = nnz(v_ac_phys  < vmin - tol); n_over_ac_tol  = nnz(v_ac_phys  > vmax + tol);

    fprintf('\n  violations (raw)      : lin under/over = %d/%d   AC under/over = %d/%d\n', ...
        n_under_lin, n_over_lin, n_under_ac, n_over_ac);
    fprintf('  violations (tol=%.0e)  : lin under/over = %d/%d   AC under/over = %d/%d\n', ...
        tol, n_under_lin_tol, n_over_lin_tol, n_under_ac_tol, n_over_ac_tol);

    ac_feasible = (n_under_ac == 0) && (n_over_ac == 0);
    fprintf('\n  NONLINEAR AC VALIDATION = %s\n', ternary(ac_feasible, 'PASSED', 'FAILED'));

    % ---- Task 6: CASE A/B/C scientific diagnosis --------------------------
    % Worst-violation magnitude vs. this point's OWN linearization error is
    % the natural yardstick for "small, boundary-sensitivity" (CASE B) vs.
    % "material" (CASE C) -- not an arbitrary absolute cutoff.
    under_viol_mag = max([0; vmin - v_ac_phys(v_ac_phys < vmin)]);
    over_viol_mag  = max([0; v_ac_phys(v_ac_phys > vmax) - vmax]);
    worst_viol_mag = max(under_viol_mag, over_viol_mag);
    max_err = max(abs(err));

    if ac_feasible
        case_class = 'A';
        fprintf('  CASE A: all nonlinear AC voltages satisfy [%.2f, %.2f] pu.\n', vmin, vmax);
        margin_reco = false;
    elseif worst_viol_mag <= max(3*max_err, 0.005)
        case_class = 'B';
        fprintf(['  CASE B: small nonlinear violation(s) near the boundary --\n' ...
            '    worst violation magnitude = %.6e pu (vs. this point''s max|error| = %.6e pu)\n' ...
            '    -> interpreted as boundary sensitivity / linearization-model mismatch,\n' ...
            '       NOT a material feasibility failure of the selected equilibrium.\n'], ...
            worst_viol_mag, max_err);
        margin_reco = true;
    else
        case_class = 'C';
        fprintf(['  CASE C: MATERIAL nonlinear violation -- worst violation magnitude = %.6e pu,\n' ...
            '    far exceeding this point''s own max|error| = %.6e pu. STOP: diagnose model\n' ...
            '    mismatch (sign/index mapping, baseline reproduction, tangent point, deviation\n' ...
            '    magnitude, stale exported interface) before drawing conclusions.\n'], ...
            worst_viol_mag, max_err);
        margin_reco = false;
    end

    replayed(end+1) = struct('label', pt.label, 'x', x, 'ac', stage.ac, ...
        'v_lin', v_lin_phys, 'v_ac', v_ac_phys, 'err', err, 'ac_feasible', ac_feasible, ...
        'case_class', case_class, 'worst_viol_mag', worst_viol_mag, ...
        'margin_reco', margin_reco); %#ok<AGROW>

    summary_rows(end+1,:) = make_summary_row(pt.label, v_lin_phys, v_ac_phys, vmin, vmax); %#ok<AGROW>
    busphase_rows = [busphase_rows; make_busphase_rows(pt.label, v_lin_phys, v_ac_phys, ...
        phys_full_idx, bus_labels, phase_labels, vmin, vmax)]; %#ok<AGROW>
end

%% ------------------------------------------------------------------------
%  Task 5: write CSV / summary outputs (never overwrite: back up first)
%  ------------------------------------------------------------------------

summary_file    = fullfile(results_dir, 'ac_validation_summary.csv');
busphase_file   = fullfile(results_dir, 'ac_validation_busphase.csv');
injection_file  = fullfile(results_dir, 'ac_validation_injections.csv');

backup_if_exists(summary_file);
backup_if_exists(busphase_file);
backup_if_exists(injection_file);

write_csv(summary_file, ...
    {'case','PF_conv','Vmin_lin','Vmin_AC','Vmax_lin','Vmax_AC','viol_lin','viol_AC','max_abs_err','rms_err'}, ...
    summary_rows);

write_csv(busphase_file, ...
    {'case','bus','phase','V_linear','V_AC','error','abs_error','under_linear','over_linear','under_AC','over_AC'}, ...
    busphase_rows);

write_csv(injection_file, ...
    {'case','bus','phase','type','delta_pu','delta_physical','DER0','DER_new'}, ...
    injection_rows);

fprintf('\nWrote:\n  %s\n  %s\n  %s\n', summary_file, busphase_file, injection_file);

%% ------------------------------------------------------------------------
%  FINAL AC VALIDATION block
%  ------------------------------------------------------------------------

get_feas = @(lbl) any(strcmp({replayed.label}, lbl)) && ...
    all(replayed(strcmp({replayed.label}, lbl)).v_ac >= vmin - 1e-9) && ...
    all(replayed(strcmp({replayed.label}, lbl)).v_ac <= vmax + 1e-9);

vgne_present = any(strcmp({replayed.label}, 'constrained vGNE'));
compromise_present = any(strcmp({replayed.label}, 'selected compromise'));
highrelief_present = any(strcmp({replayed.label}, 'high-relief selected point'));

max_err_all = 0;
for i = 1:numel(replayed)
    max_err_all = max(max_err_all, max(abs(replayed(i).err)));
end

fprintf('\n============================================================\n');
fprintf('FINAL AC VALIDATION\n');
fprintf('============================================================\n');
any_case_c = any(strcmp({replayed.case_class}, 'C'));
any_margin_reco = any([replayed.margin_reco]);

fprintf('baseline reproduced = %s\n', ternary(baseline_ok, 'YES', 'NO'));
fprintf('vGNE AC-feasible = %s\n', ternary(vgne_present, ternary(get_feas('constrained vGNE'), 'YES', 'NO'), 'NOT REPLAYED (see above)'));
fprintf('compromise AC-feasible = %s\n', ternary(compromise_present, ternary(get_feas('selected compromise'), 'YES', 'NO'), 'NOT REPLAYED -- x vector unavailable'));
fprintf('high-relief AC-feasible = %s\n', ternary(highrelief_present, ternary(get_feas('high-relief selected point'), 'YES', 'NO'), 'NOT REPLAYED -- x vector unavailable'));
fprintf('max voltage error across validated equilibria = %.6e pu\n', max_err_all);
if any_case_c
    fprintf('recommend security tightening = NO -- MATERIAL (CASE C) mismatch present, diagnose first (see above)\n');
elseif any_margin_reco
    fprintf('recommend security tightening = YES -- boundary-sensitivity (CASE B) observed; suggested Vmin_model >= %.4f pu (vmin + ~max linearization error)\n', ...
        vmin + max_err_all);
else
    fprintf('recommend security tightening = NO -- all validated points CASE A (nonlinear AC feasible as-is)\n');
end
fprintf('AC VALIDATION CLOSED = %s\n', ternary(numel(points) == 3, 'YES', 'NO -- POINTS B/C PENDING SOURCE DATA'));
fprintf('============================================================\n');

% 'points' is intentionally not saved: its x_map fields are containers.Map
% objects (input bookkeeping only); 'replayed' already holds everything
% (x as a plain vector, ac, v_lin, v_ac, err, case_class, ...) needed
% downstream, without a Map-serialization dependency.
save(fullfile(results_dir, 'ac_validation_workspace.mat'), 'replayed', ...
    'summary_rows', 'busphase_rows', 'injection_rows', 'baseline_ok', '-v7');

% Helper functions used above (make_summary_row.m, make_busphase_rows.m,
% backup_if_exists.m, write_csv.m) are standalone files at the repo root,
% following this repo's existing convention (rw.m, ternary.m, bracket.m, ...).
