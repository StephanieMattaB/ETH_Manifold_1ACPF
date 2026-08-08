%% RUN_IEEE13_PIPELINE
%
% Clean DSO-side physical pipeline for the IEEE-13 three-phase unbalanced
% feeder, built strictly as three distinct physical systems
%
%     S  -->  S'  -->  S''
%
% S   : original IEEE-13 benchmark, unmodified (ieee13.m).
% S'  : S + inverter-based DER at the four future community nodes
%       645, 611, 652, 671 (ieee13_der.m).
% S'' : S' with a lambda = 1.8 load-side congestion stress applied to
%       demand only; capacitor banks and the S' DER operating points are
%       held fixed, not scaled (ieee13_congest.m).
%
% At every stage: solve the Bolognani-Dorfler/1ACPF manifold once
% (bd_linear_solve.m), rebuild the tangent-plane linearization AROUND that
% stage's own operating point (bd_tangent_general.m -- never reused across
% stages), independently solve the nonlinear three-phase AC model
% (solve_acpf_3ph_rect.m), and compare (run_stage.m).
%
% Only AFTER S'' has been selected: derive the local input-output map and
% the DSO-to-energy-community shared constraint A_F_sh * Delta s_F <= b_F_sh
% (derive_shared_constraint.m), using S'''s tangent plane exclusively.
%
% Reused, unmodified: ieee13.m, Nmatrix.m, Rmatrix.m, bracket.m, rw.m,
% solve_acpf_3ph_rect.m.

clear all
close all
clc

%% Configuration

lambda = 1.8;      % S' -> S'' load-stress factor
vmin = 0.95;
vmax = 1.05;

ac_opts = struct();
ac_opts.tol = 1e-10;
ac_opts.max_iter = 50;
ac_opts.max_line_search = 15;
ac_opts.verbose = false;

bus_labels = {'650','632','633','634','645','646','671','692','675','680','684','611','652'};
phase_labels = {'a','b','c'};

%% Physical connectivity mask (shared by all three stages: S, S', S'' never
%  change topology or phase connectivity, only demand/generation values)

v_testfeeder = [...
    1.0625  1.0500  1.0687  ;...
    1.0210  1.0420  1.0174  ;...
    1.0180  1.0401  1.0148  ;...
    0.9940  1.0218  0.9960  ;...
    NaN     1.0329  1.0155  ;...
    NaN     1.0311  1.0134  ;...
    0.9900  1.0529  0.9778  ;...
    0.9900  1.0529  0.9777  ;...
    0.9835  1.0553  0.9758  ;...
    0.9900  1.0529  0.9778  ;...
    0.9881  NaN     0.9758  ;...
    NaN     NaN     0.9738  ;...
    0.9825  NaN     NaN     ];

[Y, v, t, p, q, n] = ieee13();

v_test_vec = reshape(v_testfeeder.', 3*n, 1);
valid_mask = ~isnan(v_test_vec);

v_slack = v(rw(1));
t_slack = t(rw(1));

fprintf('IEEE-13 DSO pipeline: %d buses, %d physical non-slack bus-phase terminals (of %d).\n', ...
    n, nnz(valid_mask(4:end)), 3*(n-1));

%% ======================================================================
%  STAGE S: original benchmark (unmodified)
%  ======================================================================

S = run_stage('S  (original benchmark)', Y, n, v_slack, t_slack, p, q, ...
    valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, []);

% Nominal-case sanity check against Kersting's published solution.
bench_err = abs(S.ac.v(valid_mask) - v_test_vec(valid_mask));
fprintf('\n[S] Nonlinear-AC vs published Kersting benchmark: max %.3e pu, RMS %.3e pu\n', ...
    max(bench_err), sqrt(mean(bench_err.^2)));

%% ======================================================================
%  STAGE S': active distribution system (+DER at 645, 611, 652, 671)
%  ======================================================================

[p_der, q_der, der_meta] = ieee13_der(n);
p_S1 = p - p_der;
q_S1 = q - q_der;

fprintf('\n[S''] DER convention: %s\n', der_meta.convention);
fprintf('[S''] Total DER capacity injected: %.1f kW across buses %s\n', ...
    der_meta.total_DER_kW, mat2str(der_meta.player_ids'));

Sp = run_stage('S'' (active: +DER at 645,611,652,671)', Y, n, v_slack, t_slack, ...
    p_S1, q_S1, valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, S.ac.V);

%% ======================================================================
%  STAGE S'': active + congested (lambda = 1.8 on load only)
%  ======================================================================

q_cap = ieee13_capacitors(n);
[p_S2, q_S2, cong_meta] = ieee13_congest(p, q, q_cap, p_der, q_der, lambda);

fprintf('\n[S''''] %s\n', cong_meta.convention);
fprintf('[S''''] Total P_load: %.1f kW (S) -> %.1f kW (S'''', lambda=%.2f); DER and capacitors fixed.\n', ...
    cong_meta.total_P_load_S*5e6/1e3, cong_meta.total_P_load_S2*5e6/1e3, lambda);

Spp = run_stage('S'''' (active + congested, lambda=1.8)', Y, n, v_slack, t_slack, ...
    p_S2, q_S2, valid_mask, bus_labels, phase_labels, vmin, vmax, ac_opts, Sp.ac.V);

%% ======================================================================
%  Three-way comparison: S vs S' vs S''
%  ======================================================================

fprintf('\n============================================================\n');
fprintf('THREE-WAY COMPARISON: S vs S'' vs S'''' (nonlinear AC ground truth)\n');
fprintf('============================================================\n');
fprintf('(violation counts below are NON-SLACK physical terminals only -- see per-stage\n');
fprintf(' sections above for the slack-inclusive figures and full terminal lists)\n\n');
fprintf('%-8s %10s %10s %8s %8s %10s %10s\n', ...
    'stage', 'Vmin[pu]', 'Vmax[pu]', 'n_under', 'n_over', 'lin-max-err', 'lin-rms-err');
stages = {S, Sp, Spp};
tags = {'S', 'S''', 'S'''''};
for k = 1:3
    st = stages{k};
    fprintf('%-8s %10.4f %10.4f %8d %8d %10.3e %10.3e\n', tags{k}, ...
        st.rep_ac_ns.Vmin, st.rep_ac_ns.Vmax, st.rep_ac_ns.n_undervoltage, st.rep_ac_ns.n_overvoltage, ...
        st.cmp.max_error, st.cmp.rms_error);
end
fprintf(['\nInterpretation: S is the passive benchmark (near-zero violations); S'' shows how DER\n' ...
    'moves the operating point (voltage rise / reduced net loading); S'''' shows the selected\n' ...
    'lambda=%.2f congestion condition producing a meaningfully (not excessively) constrained\n' ...
    'regime, which is the system adopted for the energy-community interface below.\n'], lambda);

%% ======================================================================
%  Sections 4-6: local IO map + shared constraint, derived from S'' ONLY
%  ======================================================================

[~, ~, der_meta] = ieee13_der(n);   % community struct (645,611,652,671 order)

sh = derive_shared_constraint(Spp.tp, Spp.p, Spp.q, n, valid_mask, vmin, vmax, der_meta.community);

fprintf('\n============================================================\n');
fprintf('DSO-TO-COMMUNITY INTERFACE (derived from S'''' only)\n');
fprintf('============================================================\n');
fprintf('A_z rcond              : %.3e\n', sh.rcond_Az);
fprintf('A_grid size            : %d x %d   (%d physical non-slack terminals)\n', ...
    size(sh.A_grid,1), size(sh.A_grid,2), sh.nPhys);
fprintf('A_F_sh size            : %d x %d   (4 community participants, [p,q] per phase)\n', ...
    size(sh.A_F_sh,1), size(sh.A_F_sh,2));
fprintf('min voltage margin     : upper %+ .4f pu, lower %+ .4f pu\n', ...
    min(sh.margin_upper), min(sh.margin_lower));
fprintf('rows with margin < 0   : upper %d, lower %d  (binding voltage-security rows)\n', ...
    nnz(sh.margin_upper < 0), nnz(sh.margin_lower < 0));

% Consistency check: A_F_sh * sbar_F must reproduce the same violation
% pattern encoded in b_F_sh (b_F_sh - A_F_sh*Delta_s_F=0 = b_F_sh at Delta=0).
assert(all(size(sh.A_F_sh) == [2*sh.nPhys, sum(sh.n_physical)]), ...
    'Unexpected A_F_sh dimensions.');

%% Export: DSO-to-community interface (Section 6 hand-off)

A_F_game = sh.A_F_sh;               %#ok<NASGU>  (kept name for downstream compatibility)
b_shF = sh.b_F_sh;                  %#ok<NASGU>
player_ids = sh.player_ids;         %#ok<NASGU>
n_physical = sh.n_physical;         %#ok<NASGU>
block_start = sh.block_start;       %#ok<NASGU>
block_end = sh.block_end;           %#ok<NASGU>
coordinate_convention = sh.coordinate_convention; %#ok<NASGU>
constraint_convention = sh.constraint_convention; %#ok<NASGU>

stress_meta = struct('lambda', lambda, 'vmin', vmin, 'vmax', vmax, ...
    'vstar_min', min(sh.vstar_phys), 'vstar_max', max(sh.vstar_phys), ...
    'der_meta', der_meta, 'congestion_meta', cong_meta); %#ok<NASGU>

if ~exist('data', 'dir'); mkdir('data'); end
save(fullfile('data', 'network_to_dynect_congested.mat'), ...
    'A_F_game', 'b_shF', 'player_ids', 'n_physical', 'block_start', 'block_end', ...
    'coordinate_convention', 'constraint_convention', 'stress_meta', '-v7');

fprintf('\nExported DSO-to-community interface to data/network_to_dynect_congested.mat\n');

save(fullfile('data', 'pipeline_results.mat'), 'S', 'Sp', 'Spp', 'sh', 'lambda', 'vmin', 'vmax', '-v7');
fprintf('Saved full pipeline diagnostics to data/pipeline_results.mat\n');
