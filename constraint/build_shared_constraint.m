function sh = build_shared_constraint(map, y_lo, y_hi, FL_buses, u_hat_P)
% BUILD_SHARED_CONSTRAINT  Turn the affine map y = S*u+m into linear constraints
%                          A_sh * x <= b_sh for the FL-only GNE game.
%
%   sh = build_shared_constraint(map, y_lo, y_hi, FL_buses, u_hat_P)
%
%   Implements sections 8, 10, 11 of the bridge document.
%
%   INPUTS
%     map        struct from build_io_map
%     y_lo, y_hi n_y x 1 lower/upper bounds on monitored outputs y
%     FL_buses   bus indices of strategic flexible (FL) players
%     u_hat_P    baseline (fixed) inputs of the protected (CL/NCL) nodes,
%                ordered to match the P-columns of u (see below). Optional; if
%                omitted, all non-FL inputs are held at map.u_star.
%
%   INPUT ORDERING.  Recall u = [P(nonslack); Q(nonslack)]. Each FL bus b
%   contributes TWO input entries: P_b and Q_b. A player's decision skeleton is
%   x_i = [p_i, q_i, gamma_i]; the selector Pi_i = [[1 0 0];[0 1 0]] maps
%   x_i -> (p_i,q_i) and drops gamma_i (section 11). gamma never enters physics.
%
%   OUTPUT struct `sh`:
%     A_sh_grid, b_sh_grid   grid-level  [S;-S] u <= [y_hi-m; -y_lo+m]
%     A_F, b_sh_prime        FL-block matrix and rhs after absorbing protected load
%     A_sh_blocks            cell array, one m_sh x 3 block per FL player
%     b_sh                   == b_sh_prime  (final rhs handed to solver)
%     FL_buses, m_sh
%
%   Final game constraint:   sum_i  A_sh_blocks{i} * x_i  <=  b_sh

    S = map.S;  m = map.m;
    n_y = size(S,1);
    m_sh = 2*n_y;

    % ---- section 8: grid-level two-sided -> one-sided rows -------------------
    A_sh_grid = [S; -S];
    b_sh_grid = [y_hi - m; -y_lo + m];

    % ---- input column bookkeeping -------------------------------------------
    % u = [P(nonslack); Q(nonslack)];  n_ns non-slack buses.
    nonslack = map.nonslack(:)';
    n_ns = numel(nonslack);
    colP = @(b) find(nonslack==b,1);          % column of P_b within u
    colQ = @(b) n_ns + find(nonslack==b,1);   % column of Q_b within u

    % FL input columns (p and q of every FL bus)
    FL_cols = [];
    for b = FL_buses(:)'
        cP = colP(b); cQ = colQ(b);
        if isempty(cP), error('FL bus %d is slack or invalid.', b); end
        FL_cols = [FL_cols, cP, cQ]; %#ok<AGROW>
    end
    P_cols = setdiff(1:size(S,2), FL_cols);   % protected (fixed) input columns

    S_F = S(:, FL_cols);
    S_P = S(:, P_cols);
    A_F = [S_F; -S_F];
    A_P = [S_P; -S_P];

    % ---- section 10: absorb protected load into rhs -------------------------
    if nargin < 5 || isempty(u_hat_P)
        u_hat_P = map.u_star(P_cols);         % default: baseline
    end
    b_sh_prime = b_sh_grid - A_P * u_hat_P;    % == b_sh_grid - A_P*u_hat_P

    % ---- section 11: per-player blocks with gamma column = 0 -----------------
    Pi = [1 0 0; 0 1 0];                       % (p,q) selector; drops gamma
    A_sh_blocks = cell(1, numel(FL_buses));
    for j = 1:numel(FL_buses)
        b = FL_buses(j);
        cols_j = [find(nonslack==b,1)*0 + colP(b), colQ(b)]; % [P_b, Q_b] cols
        S_Fj   = S(:, cols_j);                 % n_y x 2
        A_Fj   = [S_Fj; -S_Fj];                % m_sh x 2
        A_sh_blocks{j} = A_Fj * Pi;            % m_sh x 3, third column all zero
    end

    sh = struct();
    sh.A_sh_grid = A_sh_grid;  sh.b_sh_grid = b_sh_grid;
    sh.A_F = A_F;  sh.A_P = A_P;  sh.b_sh_prime = b_sh_prime;  sh.b_sh = b_sh_prime;
    sh.A_sh_blocks = A_sh_blocks;
    sh.FL_buses = FL_buses(:);  sh.FL_cols = FL_cols;  sh.P_cols = P_cols;
    sh.m_sh = m_sh;  sh.u_hat_P = u_hat_P;  sh.Pi = Pi;
end
