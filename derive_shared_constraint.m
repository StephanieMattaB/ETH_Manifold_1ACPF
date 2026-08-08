function sh = derive_shared_constraint(tp, p_star, q_star, n, valid_mask, vmin, vmax, community)
%DERIVE_SHARED_CONSTRAINT Sections 4-6 of the DSO pipeline: local
%input-output sensitivity map, grid-wide voltage-security constraint, and
%its projection onto the four energy-community participants.
%
%   sh = DERIVE_SHARED_CONSTRAINT(tp, p_star, q_star, n, valid_mask, ...
%                                  vmin, vmax, community)
%
%   MUST be called with the tangent plane (tp, from bd_tangent_general.m)
%   built at S'''s NONLINEAR AC-CONVERGED operating point -- NOT the
%   one-shot manifold/linear estimate. Section 3 explicitly adopts S'' "as
%   the operating system" only AFTER the manifold-vs-AC comparison; the
%   physically selected point v*'' is therefore the validated AC solution.
%   Anchoring here at the fast linear estimate instead silently changes
%   which rows the exported constraint reports as violated (see the
%   invariant check the caller is expected to run against independently
%   computed AC violation counts).
%
%   INPUTS
%   tp         : struct from bd_tangent_general.m, evaluated at S'''s
%                NONLINEAR AC operating point (v_star = tp.v_star = AC v,
%                t_star = tp.t_star = AC theta).
%   p_star,q_star : 3n x 1 net-consumption specification of S'' (the same
%                vectors passed into bd_linear_solve.m for S'').
%   valid_mask : 3n x 1 logical, physically connected terminals (S is
%                topology-invariant across S/S'/S'', so this is the same
%                mask used throughout).
%   vmin, vmax : voltage-security bounds.
%   community  : community struct array from ieee13_der.m's meta.community
%                (order 645, 611, 652, 671 -- defines A_F_sh's column order).
%
%   OUTPUT (struct sh) -- Sections 4, 5, 6 in order:
%     S_v, S_theta      : Delta z_ns = [S_v;S_theta] * Delta s_ns   (Sec 4)
%     m_v               : affine offset, v_ns = S_v*s_ns + m_v
%     A_grid, b_grid     : full grid-wide voltage constraint             (Sec 5)
%     A_F_sh, b_F_sh     : community-projected constraint, columns
%                          reordered by community participant            (Sec 6)
%     player_*          : bookkeeping for the 4 participants' column blocks

    nBus = n; nPhase = 3; nPhi = nBus*nPhase;

    idx.v = 1:nPhi; idx.theta = nPhi+(1:nPhi); idx.p = 2*nPhi+(1:nPhi); idx.q = 3*nPhi+(1:nPhi);
    idx.slackPhase = 1:3; idx.nonSlackPhase = 4:nPhi;
    idx.nonSlack_v = idx.v(idx.nonSlackPhase); idx.nonSlack_theta = idx.theta(idx.nonSlackPhase);
    nNonSlackPhi = numel(idx.nonSlackPhase);

    % --- Section 4: local input-output map, from tp only ------------------
    A_net = tp.Amat(1:2*nPhi, :);
    rowsP = 1:nPhi; rowsQ = nPhi+(1:nPhi);
    rowsNonSlack = [rowsP(idx.nonSlackPhase), rowsQ(idx.nonSlackPhase)];
    A_red = A_net(rowsNonSlack, :);

    colsZ = [idx.nonSlack_v, idx.nonSlack_theta];
    colsS = [idx.p(idx.nonSlackPhase), idx.q(idx.nonSlackPhase)];

    A_z = A_red(:, colsZ);
    A_s = A_red(:, colsS);

    assert(rank(A_z) == size(A_z,2), 'derive_shared_constraint:singularAz', 'A_z is rank deficient.');
    rcond_Az = rcond(A_z);
    if rcond_Az < 1e-8
        warning('derive_shared_constraint:illConditioned', 'A_z poorly conditioned: rcond = %.3e.', rcond_Az);
    end

    S_zs = -(A_z \ A_s);
    S_v     = S_zs(1:nNonSlackPhi, :);
    S_theta = S_zs(nNonSlackPhi+1:end, :);

    vstar_ns = tp.v_star(idx.nonSlackPhase);
    sstar_ns = [p_star(idx.nonSlackPhase); q_star(idx.nonSlackPhase)];
    m_v = vstar_ns - S_v*sstar_ns;

    valid_ns = valid_mask(idx.nonSlackPhase);
    nPhys = nnz(valid_ns);

    S_v_phys   = S_v(valid_ns, :);
    m_v_phys   = m_v(valid_ns);
    vstar_phys = vstar_ns(valid_ns);

    % --- Section 5: grid-wide voltage-security constraint ------------------
    A_grid = [S_v_phys; -S_v_phys];
    b_grid = [vmax*ones(nPhys,1) - vstar_phys; vstar_phys - vmin*ones(nPhys,1)];

    % --- Section 6: projection onto the energy community -------------------
    phase.a = 1; phase.b = 2; phase.c = 3;
    busPhaseIdx   = @(bus, ph) 3*(bus - 1) + ph;
    toNonSlackIdx = @(fullIdx) fullIdx - 3;

    nPlayers = numel(community);
    idxF_phase = cell(1, nPlayers);
    for k = 1:nPlayers
        idxF_phase{k} = toNonSlackIdx(busPhaseIdx(community(k).bus, community(k).phase));
    end
    idxF_all = [idxF_phase{:}];
    assert(numel(unique(idxF_all)) == numel(idxF_all), 'Duplicate flexible bus-phase indices.');

    idxF_p = idxF_all;
    idxF_q = nNonSlackPhi + idxF_all;
    cols_F = [idxF_p, idxF_q];

    sbar_ns = sstar_ns;
    A_F = A_grid(:, cols_F);
    sbar_F = sbar_ns(cols_F);   % kept for reporting the players' S'' baseline injection only

    % Delta s_P = 0 (rest of feeder held fixed at S''): in the PURE DEVIATION
    % formulation this drops the A_P*Delta_s_P term from the LHS and leaves
    % the RHS untouched -- per spec, b_F,sh = b_grid EXACTLY (Section 6):
    %     A_F*Delta_s_F + A_P*Delta_s_P <= b_grid,  Delta_s_P = 0
    %     => A_F*Delta_s_F <= b_grid  =>  b_F,sh := b_grid.
    % (A previous version of this code subtracted A_F*sbar_F + A_P*sbar_P
    % here, which is wrong: sbar_F/sbar_P are ABSOLUTE baseline injections,
    % not deviations, and subtracting them corrupted b_F,sh by up to ~0.19 pu
    % relative to b_grid -- caught by comparing against the AC ground truth.)
    b_F_sh = b_grid;

    % Reorder columns per-player as [p_phases, q_phases] concatenated in
    % community order (645, 611, 652, 671), derived programmatically from
    % `community` rather than hand-written column numbers.
    nF = numel(idxF_all);
    reorder = zeros(1, 2*nF);
    block_start = zeros(nPlayers,1); block_end = zeros(nPlayers,1);
    cursor = 0; pOffset = 0;
    for k = 1:nPlayers
        php = numel(idxF_phase{k});
        pLocal = (pOffset+1):(pOffset+php);
        qLocal = nF + pLocal;
        reorder(cursor+1:cursor+2*php) = [pLocal, qLocal];
        block_start(k) = cursor+1;
        block_end(k)   = cursor+2*php;
        cursor = cursor + 2*php;
        pOffset = pOffset + php;
    end

    A_F_sh = A_F(:, reorder);

    sh.idx = idx;
    sh.S_v = S_v; sh.S_theta = S_theta; sh.m_v = m_v;
    sh.rcond_Az = rcond_Az;
    sh.valid_ns = valid_ns; sh.nPhys = nPhys; sh.vstar_phys = vstar_phys;
    sh.vmin = vmin; sh.vmax = vmax;
    sh.A_grid = A_grid; sh.b_grid = b_grid;
    sh.A_F_sh = A_F_sh; sh.b_F_sh = b_F_sh;
    sh.sbar_F = sbar_F(reorder);   % players' S'' baseline injection, same column order as A_F_sh
    sh.player_ids = [community.id]';
    sh.n_physical = 2*cellfun(@numel, idxF_phase)';
    sh.block_start = block_start; sh.block_end = block_end;
    sh.margin_upper = vmax - vstar_phys;
    sh.margin_lower = vstar_phys - vmin;
    sh.coordinate_convention = 'physical variables are active/reactive DEVIATIONS from the S'' baseline';
    sh.constraint_convention = 'A_F_sh * Delta_s_F <= b_F_sh, columns ordered per player [645,611,652,671], each as [p_phases, q_phases]';
end
