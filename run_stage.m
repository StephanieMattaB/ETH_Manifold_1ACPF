function stage = run_stage(label, Y, n, v_slack, t_slack, p, q, valid_mask, ...
    bus_labels, phase_labels, vmin, vmax, ac_opts, V0_ac)
%RUN_STAGE Execute the full per-stage DSO-pipeline procedure (Sections
%1.1-1.5 / 2.1-2.5 / 3.1-3.5 of the spec, identical machinery each time):
%
%   1-2. one-shot Bolognani-Dorfler manifold solve         -> x*
%   3.   fresh tangent-plane linearization AROUND x*        (never reused)
%   4.   independent nonlinear three-phase AC solve
%   5.   manifold-vs-AC comparison
%   + full voltage-security diagnostics for both solutions.
%
%   stage = RUN_STAGE(label, Y, n, v_slack, t_slack, p, q, valid_mask, ...
%                      bus_labels, phase_labels, vmin, vmax, ac_opts, V0_ac)
%
%   p, q       : 3n x 1 net-consumption specification for THIS stage.
%   V0_ac      : optional 3n x 1 complex initial guess for the AC solver
%                (pass [] for the flat/slack-repeated default).

    fprintf('\n============================================================\n');
    fprintf('STAGE %s\n', label);
    fprintf('============================================================\n');

    % --- 1-2: manifold solve -------------------------------------------------
    [v_lin, t_lin] = bd_linear_solve(Y, n, p, q, v_slack, t_slack);

    % --- 3: fresh tangent-plane linearization around x* ----------------------
    tp = bd_tangent_general(Y, n, v_lin, t_lin);

    % Manifold residual: how far the (v_lin,t_lin) point sits from the TRUE
    % nonlinear manifold, i.e. || -Sstar(v_lin,t_lin) - specified S ||  at the
    % non-slack physical rows (0 would mean x* is exactly on-manifold).
    S_manifold_spec = -(p + 1j*q);
    nonSlackFull = rw(2:n);
    tp.manifold_residual = norm(tp.Sstar(nonSlackFull) - S_manifold_spec(nonSlackFull), inf);

    % --- 4: independent nonlinear AC solve ------------------------------------
    Vslack = v_slack .* exp(1j*t_slack);
    if nargin < 14 || isempty(V0_ac)
        V0_ac = repmat(Vslack, n, 1);
    end
    Sspec = -(p + 1j*q);
    Sspec(1:3) = 0;   % slack rows unused by solve_acpf_3ph_rect

    ac = solve_acpf_3ph_rect(Y, (1:3)', Vslack, Sspec, V0_ac, ac_opts);
    if ~ac.converged
        error('run_stage:acNotConverged', ...
            'Stage %s: nonlinear AC power flow did not converge (max mismatch %.3e).', ...
            label, ac.max_mismatch);
    end

    % --- 5: voltage-security diagnostics + manifold-vs-AC comparison ---------
    % Slack-bus phases (1:3) are a fixed regulator SETTING (an input, not a
    % solved network state); report them separately from genuine non-slack
    % violations so a chronically-set slack magnitude (e.g. 1.0625/1.0687 pu
    % here) is never confused with a network-security issue.
    slack_mask = false(size(valid_mask)); slack_mask(1:3) = true;
    valid_mask_ns = valid_mask & ~slack_mask;

    rep_lin = voltage_report(v_lin, valid_mask, bus_labels, phase_labels, vmin, vmax);
    rep_ac  = voltage_report(ac.v,  valid_mask, bus_labels, phase_labels, vmin, vmax);
    rep_lin_ns = voltage_report(v_lin, valid_mask_ns, bus_labels, phase_labels, vmin, vmax);
    rep_ac_ns  = voltage_report(ac.v,  valid_mask_ns, bus_labels, phase_labels, vmin, vmax);
    cmp     = compare_voltages(v_lin, ac.v, valid_mask, bus_labels, phase_labels);

    fprintf('\n-- Nonlinear AC --\n');
    fprintf('  converged            : %d (%d Newton iterations, max mismatch %.3e)\n', ...
        ac.converged, ac.iterations, ac.max_mismatch);
    fprintf('  Vmin                 : %.6f pu  at %s.%s\n', rep_ac.Vmin, rep_ac.Vmin_bus, rep_ac.Vmin_phase);
    fprintf('  Vmax                 : %.6f pu  at %s.%s\n', rep_ac.Vmax, rep_ac.Vmax_bus, rep_ac.Vmax_phase);
    fprintf('  under/over/total viol: %d / %d / %d  (of %d physical terminals, incl. slack)\n', ...
        rep_ac.n_undervoltage, rep_ac.n_overvoltage, rep_ac.n_violations, rep_ac.n_physical);
    fprintf('  ...of which non-slack : %d / %d / %d  (of %d non-slack physical terminals)\n', ...
        rep_ac_ns.n_undervoltage, rep_ac_ns.n_overvoltage, rep_ac_ns.n_violations, rep_ac_ns.n_physical);
    if rep_ac.n_violations > 0
        fprintf('  violated terminals    : %s\n', strjoin(rep_ac.violated_terminals, ', '));
    end

    fprintf('\n-- Manifold / linear (1ACPF one-shot) --\n');
    fprintf('  Vmin                 : %.6f pu  at %s.%s\n', rep_lin.Vmin, rep_lin.Vmin_bus, rep_lin.Vmin_phase);
    fprintf('  Vmax                 : %.6f pu  at %s.%s\n', rep_lin.Vmax, rep_lin.Vmax_bus, rep_lin.Vmax_phase);
    fprintf('  under/over/total viol: %d / %d / %d  (of %d physical terminals, incl. slack)\n', ...
        rep_lin.n_undervoltage, rep_lin.n_overvoltage, rep_lin.n_violations, rep_lin.n_physical);
    fprintf('  ...of which non-slack : %d / %d / %d  (of %d non-slack physical terminals)\n', ...
        rep_lin_ns.n_undervoltage, rep_lin_ns.n_overvoltage, rep_lin_ns.n_violations, rep_lin_ns.n_physical);
    fprintf('  manifold residual     : %.3e pu (||spec - self-consistent S|| at x*)\n', tp.manifold_residual);
    fprintf('  tangent-plane rcond   : %.3e (full assembled system)\n', tp.rcond_top);

    fprintf('\n-- Manifold vs nonlinear AC --\n');
    fprintf('  max |v_lin - v_AC|    : %.6e pu  (worst at %s.%s: %.4f vs %.4f)\n', ...
        cmp.max_error, cmp.worst_bus, cmp.worst_phase, cmp.worst_v_lin, cmp.worst_v_ac);
    fprintf('  RMS |v_lin - v_AC|    : %.6e pu\n', cmp.rms_error);

    stage.label = label;
    stage.p = p; stage.q = q;
    stage.v_lin = v_lin; stage.t_lin = t_lin; stage.tp = tp;
    stage.ac = ac;
    stage.rep_lin = rep_lin; stage.rep_ac = rep_ac; stage.cmp = cmp;
    stage.rep_lin_ns = rep_lin_ns; stage.rep_ac_ns = rep_ac_ns;
end
