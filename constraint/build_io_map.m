function map = build_io_map(Ybus, V0, A0, slack, mon_buses)
% BUILD_IO_MAP  Explicit affine input-output map  y = S*u_abs + m  from the
%               Bolognani-Dorfler tangent manifold (v1: PQ-for-everyone).
%
%   map = build_io_map(Ybus, V0, A0, slack, mon_buses)
%
%   Implements sections 5-7 of the "ACPF Tangent Manifold -> Shared Constraint"
%   bridge. v1 scope (section 9/13): the entire non-slack network is PQ-typed,
%   i.e. every non-slack (p_i, q_i) is a controllable INPUT and every non-slack
%   (Vmag_i, theta_i) is a solved UNKNOWN. Monitored outputs y are the voltage
%   MAGNITUDES at the buses listed in mon_buses.
%
%   INPUTS
%     Ybus       n x n complex bus admittance matrix
%     V0, A0     n x 1 linearization point (magnitudes, angles [rad]).
%                Pass the TRUE AC operating point for a locally exact map.
%     slack      index of the slack (REF) bus
%     mon_buses  vector of bus indices whose voltage magnitude is monitored
%
%   OUTPUT struct `map` with fields:
%     S          n_y x n_u sensitivity  d(y)/d(u)
%     m          n_y x 1  affine offset, y = S*u_abs + m
%     u_star     n_u x 1  baseline input (absolute [P_K; Q_K] at lin. point)
%     y_star     n_y x 1  baseline output (Vmag at mon_buses)
%     A_U, A_K   the (2n-2)x(2n-2) partition blocks (kept for debugging)
%     C_U, C_K   output selectors
%     u_idx      struct describing input ordering:  u = [P(nonslack); Q(nonslack)]
%     out_idx    monitored bus list
%     A_vd       full 2n x 2n tangent matrix
%
%   The input vector u is ordered   u = [ P_i ; Q_i ]  for i in nonslack buses,
%   nonslack buses taken in ascending index order.

    n = numel(V0);
    A_vd = build_Avd(Ybus, V0, A0);          % 2n x 2n, maps [dV;dTh] -> [dP;dQ]

    % ---- index sets (state layout: 1..n = V/P block, n+1..2n = Th/Q block) ----
    nonslack = setdiff(1:n, slack);
    % voltage-side UNKNOWN columns: dVmag and dTheta of non-slack buses
    U_cols = [nonslack, n + nonslack];
    % power-side INPUT rows: dP and dQ of non-slack buses
    K_rows = [nonslack, n + nonslack];

    % ---- partition:  A_vd * [dV;dTh] = [dP;dQ]  ------------------------------
    % Slack voltage cols drop (dVmag_slk = dTheta_slk = 0).
    % Slack power rows drop (they define the slack injection, an output).
    % Result:   A_U * z_U = u        with  A_U = A_vd(K_rows, U_cols)
    % In the doc's general form  A_U*xiU + A_K*xiK = 0 with A_K = -I (power = -I*s),
    % so T = -A_U\A_K = A_U\I.
    A_U = A_vd(K_rows, U_cols);              % (2n-2) x (2n-2)
    A_K = -eye(numel(K_rows));               % power enters as -I

    rc = rcond(A_U);
    if ~isfinite(rc) || rc < 1e-12
        warning('build_io_map:illcond', ...
            'A_U is ill-conditioned (rcond = %.3e). Map may be unreliable.', rc);
    end

    % ---- output selectors:  y = Vmag at mon_buses  --------------------------
    % z_U is ordered as U_cols; within it, dVmag of bus b sits at the position
    % of b inside `nonslack`. Monitored buses must be non-slack for v1.
    n_y = numel(mon_buses);
    n_u = numel(K_rows);
    C_U = zeros(n_y, n_u);
    for r = 1:n_y
        b = mon_buses(r);
        pos = find(nonslack == b, 1);
        if isempty(pos)
            error('build_io_map:monSlack', ...
                'Monitored bus %d is the slack (or invalid); pick a non-slack bus.', b);
        end
        C_U(r, pos) = 1;                    % dVmag lives in first half of z_U
    end
    C_K = zeros(n_y, n_u);                  % Vmag does not depend on u except via z_U

    % ---- the map:  S = C_K - C_U*(A_U \ A_K) ; m = y_star - S*u_star ---------
    S = C_K - C_U * (A_U \ A_K);            % linear solve, NEVER inv(A_U)

    % baseline (absolute) input and output at the linearization point
    Ucplx  = V0 .* exp(1j*A0);
    Scplx  = Ucplx .* conj(Ybus * Ucplx);  % complex power injections
    Pinj   = real(Scplx);
    Qinj   = imag(Scplx);
    u_star = [Pinj(nonslack); Qinj(nonslack)];
    y_star = V0(mon_buses);

    m = y_star - S * u_star;

    map = struct();
    map.S = S;  map.m = m;
    map.u_star = u_star;  map.y_star = y_star;
    map.A_U = A_U;  map.A_K = A_K;  map.C_U = C_U;  map.C_K = C_K;
    map.A_vd = A_vd;
    map.slack = slack;  map.nonslack = nonslack;
    map.mon_buses = mon_buses(:);
    map.u_idx = struct('order', 'u = [P(nonslack); Q(nonslack)]', ...
                       'nonslack', nonslack, 'n', n);
    map.rcond_AU = rc;
end
