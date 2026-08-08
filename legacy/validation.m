% validation.m — linearization vs nonlinear-AC comparison
%
% Prerequisites (run in order before this script):
%   1. example_ieee13.m  → puts vstar_phys, valid_ns, idx, sbar_F in workspace
%   2. run_ieee13_nonlinear_ac.m  → writes nonlinear_ac_validation_results.mat
%
% Outputs:
%   • max |v_lin − v_AC|, RMS error, worst phase
%   • AC stressed-case summary (Vmin, Vmax, violations)
%   • Per-player baseline injections (linearized convention, net consumption)

load('nonlinear_ac_validation_results.mat')   % → single struct: results

% ── Lin vs AC voltage error ───────────────────────────────────────────────────
% Map the 29 physical rows back to the 39-slot AC vector.
% idx.nonSlackPhase (4:39, 36 elements) filtered by valid_ns (29-of-36 logical).
phys_full_idx = idx.nonSlackPhase(valid_ns);      % 39-slot indices, 29 elements
v_AC_ns  = results.stressed.v(phys_full_idx);     % 29×1 AC voltages

err      = abs(vstar_phys - v_AC_ns);
[maxerr, widx] = max(err);

fprintf('=== Linearization error (lin vs nonlinear AC at %.2fx) ===\n', ...
        results.settings.stress_factor);
fprintf('  max |v_lin − v_AC| = %.6e pu\n', maxerr);
fprintf('  RMS(v_lin − v_AC)  = %.6e pu\n', sqrt(mean(err.^2)));
fprintf('  worst row %2d : v_lin = %.4f  v_AC = %.4f  Δ = %+.4f\n', ...
        widx, vstar_phys(widx), v_AC_ns(widx), vstar_phys(widx) - v_AC_ns(widx));

% ── AC stressed summary ───────────────────────────────────────────────────────
s = results.stressed.summary;
fprintf('\n=== AC stressed (%.2fx) voltage summary ===\n', ...
        results.settings.stress_factor);
fprintf('  Vmin = %.6f pu  at bus %s phase %s\n', s.Vmin, s.Vmin_bus, s.Vmin_phase);
fprintf('  Vmax = %.6f pu  at bus %s phase %s\n', s.Vmax, s.Vmax_bus, s.Vmax_phase);
fprintf('  Under-voltage: %d   Over-voltage: %d\n', ...
        s.n_undervoltage, s.n_overvoltage);

% ── Linearized summary for comparison ────────────────────────────────────────
fprintf('\n=== Linearized (example_ieee13) voltage summary at %.2fx ===\n', ...
        results.settings.stress_factor);
fprintf('  Vmin = %.6f pu\n', min(vstar_phys));
fprintf('  Vmax = %.6f pu\n', max(vstar_phys));
[~, lo_r] = min(vstar_phys);
[~, hi_r] = max(vstar_phys);
fprintf('  min at row %d, max at row %d  (row = position in 29-element phys vector)\n', ...
        lo_r, hi_r);

% ── Per-player baseline injections (linearized, net-consumption convention) ───
%
% sbar_F layout  (14 elements):
%   p-block (elements 1..7): p645_b  p645_c  p611_c  p652_a  p671_a  p671_b  p671_c
%   q-block (elements 8..14): q645_b  q645_c  q611_c  q652_a  q671_a  q671_b  q671_c
%
% Positive value = net CONSUMPTION (load convention from ieee13.m).
% Game deviation x_i = Δs_F (same convention); players reduce consumption to help voltage.

player_ids = [645,    611,   652,   671];
p_idx      = {[1 2],  [3],   [4],   [5 6 7]};
q_idx      = {[8 9],  [10],  [11],  [12 13 14]};
ph_labels  = {{'b','c'}, {'c'}, {'a'}, {'a','b','c'}};

fprintf('\n=== Per-player baseline at stressed point (linearized, net-consumption) ===\n');
for k = 1:4
    p_bar = sbar_F(p_idx{k});
    q_bar = sbar_F(q_idx{k});
    phs   = ph_labels{k};
    fprintf('  player %d:\n', player_ids(k));
    for ph = 1:numel(phs)
        fprintf('    phase %s:  p̄ = %+.5f pu   q̄ = %+.5f pu\n', ...
                phs{ph}, p_bar(ph), q_bar(ph));
    end
end
