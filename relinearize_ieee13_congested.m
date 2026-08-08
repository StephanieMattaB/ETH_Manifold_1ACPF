%% RELINEARIZE_IEEE13_CONGESTED
%
%   Standalone, additive script. Does NOT modify and does NOT call
%   example_ieee13.m or run_ieee13_nonlinear_ac.m -- it reproduces the
%   small amount of setup it needs directly, using only the existing
%   helper files (ieee13.m, Nmatrix.m, Rmatrix.m, bracket.m, rw.m), so
%   nothing in the existing pipeline is at risk of being broken.
%
%   PURPOSE
%   -------
%   example_ieee13.m builds its linear power-flow model (Amat) by
%   linearizing (implicitly, a la Bolognani-Dorfler) around the FLAT /
%   NOMINAL reference voltage profile: every bus-phase is referenced to
%   the same ideal triple aaa = [1; a; a^2] (magnitude 1, balanced
%   120-degree angles), via
%
%       RRR = bracket(kron(eye(n), diag(aaa)))
%           = bracket(diag(kron(ones(n,1), aaa)))
%
%   That is why S_v (the voltage sensitivity to injections) in that
%   script is load-INDEPENDENT: the reference point never moves, only
%   the RHS (p,q) does. This script relinearizes around the ACTUAL
%   stressed / congested operating point instead, by substituting the
%   real solved voltage phasor for the reference:
%
%       RRR = bracket(diag(Vref))
%
%   and iterating Vref -> solve -> new Vref to a fixed point (a Picard
%   iteration toward the true nonlinear AC solution), using nothing but
%   this same linear machinery -- no nonlinear solver dependency.
%
%   OUTPUT
%   ------
%   data/network_to_dynect_congested.mat, in the same format/variable
%   names as example_ieee13.m's network_to_dynect.mat, so the DyNECT /
%   OSQP side can consume either without changes. The original file is
%   never written by this script.

%% Network + baseline setup (mirrors example_ieee13.m, unmodified there)

Vbase = 4160/sqrt(3);
Sbase = 5e6;
Zbase = Vbase^2/Sbase;

[Y,v,t,p,q,n] = ieee13();

e0  = [1; zeros(n-1,1)];
a   = exp(-1j*2*pi/3);
aaa = [1; a; a^2];

VTV = [kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n), zeros(3, 3*n)];
VTT = [zeros(3, 3*n), kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n)];
PQP = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1)), zeros(3*(n-1),3*n)];
PQQ = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1))];

NNN = Nmatrix(6*n);
LLL = bracket(Y);

%% Congestion stress (must match example_ieee13.m's operating point exactly)

load_scale = 1.8;

bus_extra = containers.Map('KeyType','double','ValueType','double');

q_cap = zeros(size(q));
q_cap(rw(9))    = [200; 200; 200]*1e3/Sbase;
q_cap(rw(12,3)) = 100e3/Sbase;

p = load_scale*p;
q = load_scale*(q + q_cap) - q_cap;

for kk = keys(bus_extra)
    bb = kk{1};
    p(rw(bb)) = bus_extra(bb)*p(rw(bb));
    q(rw(bb)) = bus_extra(bb)*(q(rw(bb)) + q_cap(rw(bb))) ...
              - q_cap(rw(bb));
end

fprintf('Congestion stress: global load_scale = %.3f\n', load_scale);

Bmat = [zeros(3*n,1); zeros(3*n,1); v(rw(1)); t(rw(1)); p(rw(2:n)); q(rw(2:n))];

%% Flat/no-load-anchored baseline (same model as example_ieee13.m, recomputed
%  here only as the fixed-point iteration's starting point and as a
%  reference for the comparison printout at the end -- NOT exported).

RRR_flat  = bracket(kron(eye(n), diag(aaa)));
Amat_flat = [NNN*inv(RRR_flat)*LLL*RRR_flat, eye(6*n); VTV; VTT; PQP; PQQ];
x_flat    = Amat_flat \ Bmat;

%% ============================================================
%  Relinearization around the congested (stressed) operating point
%  ------------------------------------------------------------
%  IMPORTANT CORRECTNESS NOTES (found by checking against the source
%  paper, Bolognani & Dorfler, Allerton 2015, Proposition 1 -- both
%  bugs were caught by actually running this script, not by inspection):
%
%  (1) Amat_flat above is NOT the general implicit-linearization
%  formula -- it is the FLAT-ONLY special case (their eq. 7 / the
%  three-phase variant in Sec. III-E), valid only because at v*=1
%  (every magnitude exactly 1) the identity <diag(u*)> = R(u*) holds,
%  letting RRR = bracket(diag(u*)) stand in for R(u*) inside a
%  shunt-free simplified expression. That identity is FALSE once
%  magnitudes deviate from 1 -- exactly what happens at a real
%  (stressed) operating point. (A first attempt substituted the actual
%  operating point directly into that flat-only formula; it diverged
%  numerically, confirming the identity does not generalize.)
%
%  The GENERAL tangent-plane formula (their eq. 4-5), valid at ANY
%  point x* = (v*,theta*,p*,q*) on the power-flow manifold, is
%
%     A_x* = [ M(u*) * R(u*),  I ]     (I, not -I, because this file's
%                                       p,q are NET CONSUMPTION =
%                                       -injection; see
%                                       run_ieee13_nonlinear_ac.m)
%     M(u*) = bracket(diag(conj(Y*u*))) + bracket(diag(u*))*NNN*bracket(Y)
%     R(u*) = Rmatrix(v*, theta*)      <- their own helper, unused
%                                          (dead code) in example_ieee13.m
%
%  Verified below to reduce to Amat_flat's top-left block at the exact
%  flat point (v*=1 everywhere), to machine precision.
%
%  (2) The tangent-plane equation is A_x*(x - x*) = 0, i.e.
%  A_x* x = A_x* x*, NOT "A_x* x = 0". A_x* x* is generally NONZERO
%  away from the flat point (it collapses to 0 there only because
%  y^sh=0 for this feeder and the flat-point self-consistent injection
%  happens to be exactly 0). A first attempt kept the RHS at 0 for
%  every iterate (copying Amat_flat's Bmat unchanged) and diverged --
%  the fix is to recompute the true RHS, A_x* * x*, at each iterate
%  from the CURRENT (v*,theta*) and its self-consistent injection
%  S* = u*.*conj(Y*u*) (the exact nonlinear relation, so x* trivially
%  sits on the manifold by construction). With that fix this is
%  literally Newton's method for the power flow, expressed in this
%  paper's coordinates instead of raw P/Q-mismatch coordinates --
%  using only Rmatrix/bracket/Nmatrix/Y, nothing from
%  solve_acpf_3ph_rect.m. It converges quadratically once correct.
%
%  Convergence is checked two ways, both self-contained:
%    (1) the voltage estimate stops moving (||d[v;theta]|| small);
%    (2) the TRUE nonlinear AC mismatch V.*conj(Y*V) - Sspec (built
%        directly from Y) is small.
%% ============================================================

relin.tol_V    = 1e-9;
relin.tol_ac   = 1e-9;
relin.max_iter = 30;

nonSlackFull   = rw(2:n);                                   % full 39-slot non-slack indices
Sspec_nonslack = -(p(nonSlackFull) + 1j*q(nonSlackFull));    % injection convention

% One-time formula self-check at the EXACT flat point (v=1 everywhere,
% matching Amat_flat's own linearization point precisely) -- confirms
% M(u*)*R(u*) reduces to Amat_flat's top-left block before trusting it
% anywhere else.
v_check = ones(3*n,1);
t_check = kron(ones(n,1), angle(aaa));
Vcheck  = v_check .* exp(1j*t_check);
Mcheck  = bracket(diag(conj(Y*Vcheck))) + bracket(diag(Vcheck))*NNN*bracket(Y);
Rcheck  = Rmatrix(v_check, t_check);
selfcheck = norm(Mcheck*Rcheck - NNN*inv(RRR_flat)*LLL*RRR_flat, 'fro') / ...
    norm(NNN*inv(RRR_flat)*LLL*RRR_flat, 'fro');
fprintf('\nFormula self-check at exact flat point: relative diff = %.3e\n', selfcheck);
if selfcheck > 1e-8
    error(['General-formula construction does NOT reduce to Amat_flat''s top-left ', ...
        'block at the flat point (relative diff %.3e) -- refusing to trust the ', ...
        'congested-point relinearization below.'], selfcheck);
end

v_est = ones(3*n,1);
t_est = kron(ones(n,1), angle(aaa));
v_est(1:3) = v(rw(1));
t_est(1:3) = t(rw(1));

x_cong = x_flat;                 % fallback if the loop body never runs

fprintf('\n=== Relinearizing around the congested operating point (Newton, general formula) ===\n');

for it = 0:relin.max_iter

    Vk = v_est .* exp(1j*t_est);
    Sk = Vk .* conj(Y*Vk);                     % self-consistent injection at the current iterate
    s_cons_k = -[real(Sk); imag(Sk)];          % [p*;q*] in this file's consumption convention

    Mk = bracket(diag(conj(Y*Vk))) + bracket(diag(Vk))*NNN*bracket(Y);
    Rk = Rmatrix(v_est, t_est);
    Atop = Mk*Rk;

    Amat_cong = [Atop, eye(6*n); VTV; VTT; PQP; PQQ];

    rhs_top   = Atop*[v_est; t_est] + s_cons_k;   % = A_x* * x*  (generally nonzero away from flat)
    Bmat_cong = [rhs_top; Bmat(6*n+1:end)];       % bottom (bus-model) rows are x*-independent, unchanged

    x_cong = Amat_cong \ Bmat_cong;

    v_new = x_cong(1:3*n);
    t_new = x_cong(3*n+1:2*3*n);

    dV = norm([v_new;t_new] - [v_est;t_est], inf);

    V_new       = v_new .* exp(1j*t_new);
    Scalc_it    = V_new .* conj(Y*V_new);
    ac_mismatch = norm(Scalc_it(nonSlackFull) - Sspec_nonslack, inf);

    fprintf('  iter %2d: ||d[v;theta]||_inf = %.3e   true AC mismatch = %.3e\n', ...
        it, dV, ac_mismatch);

    v_est = v_new;
    t_est = t_new;
    relin.iterations = it;

    if dV < relin.tol_V && ac_mismatch < relin.tol_ac
        fprintf('  Converged in %d iterations.\n', it);
        relin.converged = true;
        break;
    end

    relin.converged = false;
    if it == relin.max_iter
        warning(['Relinearization Newton iteration did not converge in %d ', ...
            'iterations (||d[v;theta]||=%.3e, AC mismatch=%.3e).'], ...
            it, dV, ac_mismatch);
    end
end

relin.final_dV          = dV;
relin.final_ac_mismatch = ac_mismatch;

%% Benchmark table (only used for the valid_ns physical-row mask)

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

v_testfeeder = reshape(v_testfeeder.', 3*n, 1);

%% Run the shared sensitivity/aggregation pipeline on both models

fit_flat = build_shared_constraint(Amat_flat, Bmat, x_flat, n, v_testfeeder, ...
    load_scale, 'flat/no-load anchored (reference only, not exported)');

fit_cong = build_shared_constraint(Amat_cong, Bmat_cong, x_cong, n, v_testfeeder, ...
    load_scale, 'congested-point anchored');

%% Export only the congested-point model

export_shared_constraint(fit_cong, "network_to_dynect_congested.mat", ...
    load_scale, relin);

%% Comparison printout

fprintf('\n=== Flat (no-load) vs congested-point comparison ===\n');
fprintf('S_v range        : flat [%.4f, %.4f]   congested [%.4f, %.4f]\n', ...
    min(fit_flat.S_v(:)), max(fit_flat.S_v(:)), min(fit_cong.S_v(:)), max(fit_cong.S_v(:)));
fprintf('m_v range        : flat [%.4f, %.4f]   congested [%.4f, %.4f]\n', ...
    min(fit_flat.m_v), max(fit_flat.m_v), min(fit_cong.m_v), max(fit_cong.m_v));
fprintf('vstar_phys range : flat [%.4f, %.4f]   congested [%.4f, %.4f]\n', ...
    min(fit_flat.vstar_phys), max(fit_flat.vstar_phys), min(fit_cong.vstar_phys), max(fit_cong.vstar_phys));
fprintf('min margin_lower : flat %+.4f   congested %+.4f\n', ...
    min(fit_flat.margin_lower), min(fit_cong.margin_lower));
fprintf('min margin_upper : flat %+.4f   congested %+.4f\n', ...
    min(fit_flat.margin_upper), min(fit_cong.margin_upper));


%% ============================================================
%  Local functions (must follow all script commands in a MATLAB script)
%% ============================================================

function out = build_shared_constraint(A, b, xstar, n, v_testfeeder, load_scale, label)
%BUILD_SHARED_CONSTRAINT Voltage-sensitivity constraint + flexible-community
%aggregation, exactly as in example_ieee13.m, parametrized so it can be run
%on any linearization (A,b,xstar) -- flat/no-load or congested-point.

    nBus   = 13;
    nPhase = 3;
    nPhi   = nBus*nPhase;       % 39 bus-phase slots

    assert(isequal(size(A), [4*nPhi, 4*nPhi]), 'A must be 156 x 156.');
    assert(isequal(size(b), [4*nPhi, 1]), 'b must be 156 x 1.');
    assert(isequal(size(xstar), [4*nPhi, 1]), 'xstar must be 156 x 1.');

    % State blocks: x = [v; theta; p; q]
    idx.v     = 1:nPhi;
    idx.theta = nPhi + (1:nPhi);
    idx.p     = 2*nPhi + (1:nPhi);
    idx.q     = 3*nPhi + (1:nPhi);

    idx.slackPhase    = 1:3;
    idx.nonSlackPhase = 4:nPhi;

    idx.slack_v     = idx.v(idx.slackPhase);
    idx.slack_theta = idx.theta(idx.slackPhase);

    idx.nonSlack_v     = idx.v(idx.nonSlackPhase);
    idx.nonSlack_theta = idx.theta(idx.nonSlackPhase);
    idx.nonSlack_p      = idx.p(idx.nonSlackPhase);
    idx.nonSlack_q      = idx.q(idx.nonSlackPhase);

    nNonSlackPhi = numel(idx.nonSlackPhase);   % 36

    %% Network manifold (ROWS): 1:39 active power, 40:78 reactive power
    A_net = A(1:2*nPhi, :);
    b_net = b(1:2*nPhi);

    rowsP = 1:nPhi;
    rowsQ = nPhi + (1:nPhi);
    rowsNonSlack = [rowsP(idx.nonSlackPhase), rowsQ(idx.nonSlackPhase)];

    A_red = A_net(rowsNonSlack, :);
    b_red = b_net(rowsNonSlack);

    %% Partition reduced variables
    colsZ  = [idx.nonSlack_v, idx.nonSlack_theta];   % z_ns = [v_ns; theta_ns]
    colsS  = [idx.nonSlack_p, idx.nonSlack_q];        % s_ns = [p_ns; q_ns]
    colsZ0 = [idx.slack_v, idx.slack_theta];          % z_0  = [v_0; theta_0]

    A_z = A_red(:, colsZ);
    A_s = A_red(:, colsS);
    A_0 = A_red(:, colsZ0);

    assert(isequal(size(A_z), [2*nNonSlackPhi, 2*nNonSlackPhi]), ...
        'Unexpected dimensions for A_z.');
    assert(rank(A_z) == size(A_z,2), 'A_z is rank deficient.');
    if rcond(A_z) < 1e-8
        warning('[%s] A_z is poorly conditioned: rcond = %.3e.', label, rcond(A_z));
    end

    %% Reference operating point
    zstar_ns = [xstar(idx.nonSlack_v); xstar(idx.nonSlack_theta)];
    sstar_ns = [xstar(idx.nonSlack_p); xstar(idx.nonSlack_q)];
    zstar_0  = [xstar(idx.slack_v); xstar(idx.slack_theta)];
    vstar_ns = xstar(idx.nonSlack_v);

    %% Sensitivity map: Delta z_ns = S_zs * Delta s_ns
    S_zs = -(A_z \ A_s);
    S_v     = S_zs(1:nNonSlackPhi, :);
    S_theta = S_zs(nNonSlackPhi+1:end, :);

    %% Affine voltage map: v_ns = S_v*s_ns + m_v
    m_v = vstar_ns - S_v*sstar_ns;

    valid_ns = ~isnan(v_testfeeder(idx.nonSlackPhase));
    nPhys = nnz(valid_ns);

    S_v_phys   = S_v(valid_ns,:);
    m_v_phys   = m_v(valid_ns);
    vstar_phys = vstar_ns(valid_ns);

    fprintf('\n[%s] m_v range: [%.4f, %.4f] p.u.\n', label, min(m_v), max(m_v));

    %% Voltage constraint
    vmin = 0.95*ones(nPhys,1);
    vmax = 1.05*ones(nPhys,1);

    A_grid = [S_v_phys; -S_v_phys];
    b_grid = [vmax - m_v_phys; -vmin + m_v_phys];

    %% Validation
    fullResidual    = norm(A*xstar - b, 2);
    reducedResidual = norm(A_z*zstar_ns + A_s*sstar_ns + A_0*zstar_0 - b_red, 2);
    mapResidual     = norm(vstar_ns - (S_v*sstar_ns + m_v), inf);
    constraintResidual = A_grid*sstar_ns - b_grid;

    fprintf('[%s] Full manifold residual:    %.3e\n', label, fullResidual);
    fprintf('[%s] Reduced manifold residual: %.3e\n', label, reducedResidual);
    fprintf('[%s] Affine-map residual:       %.3e\n', label, mapResidual);
    fprintf('[%s] Max constraint violation:  %.3e\n', label, max(constraintResidual));

    %% Flexible-community bus-phase locations
    % 645 -> bus 5, phases b,c | 611 -> bus 12, phase c
    % 652 -> bus 13, phase a   | 671 -> bus 7, phases a,b,c
    phase.a = 1; phase.b = 2; phase.c = 3;
    busPhaseIdx   = @(bus, ph) 3*(bus - 1) + ph;
    toNonSlackIdx = @(fullIdx) fullIdx - 3;

    idxF_645 = toNonSlackIdx([busPhaseIdx(5, phase.b), busPhaseIdx(5, phase.c)]);
    idxF_611 = toNonSlackIdx(busPhaseIdx(12, phase.c));
    idxF_652 = toNonSlackIdx(busPhaseIdx(13, phase.a));
    idxF_671 = toNonSlackIdx([busPhaseIdx(7, phase.a), busPhaseIdx(7, phase.b), busPhaseIdx(7, phase.c)]);

    idxF_phase = [idxF_645, idxF_611, idxF_652, idxF_671];
    assert(numel(unique(idxF_phase)) == numel(idxF_phase), ...
        'Duplicate flexible bus-phase indices detected.');

    idxF_p = idxF_phase;
    idxF_q = nNonSlackPhi + idxF_phase;
    cols_F = [idxF_p, idxF_q];
    cols_all = 1:(2*nNonSlackPhi);
    cols_P = setdiff(cols_all, cols_F, 'stable');

    sbar_ns = sstar_ns;
    sbar_F = sbar_ns(cols_F);
    sbar_P = sbar_ns(cols_P);

    A_F = A_grid(:, cols_F);
    A_P = A_grid(:, cols_P);

    b_shF = b_grid - A_F*sbar_F - A_P*sbar_P;
    b_shF_check = b_grid - A_grid*sbar_ns;
    assert(norm(b_shF - b_shF_check, inf) < 1e-10, ...
        'Flexible/protected partition is inconsistent.');

    fprintf('[%s] A_F size: %d x %d   b_shF size: %d x %d\n', ...
        label, size(A_F,1), size(A_F,2), size(b_shF,1), size(b_shF,2));

    % Physical flexible-variable ordering in A_F:
    % [p645_b,p645_c,p611_c,p652_a,p671_a,p671_b,p671_c, q645_b,q645_c,q611_c,q652_a,q671_a,q671_b,q671_c]
    colsPlayer_645 = [1 2 8 9];
    colsPlayer_611 = [3 10];
    colsPlayer_652 = [4 11];
    colsPlayer_671 = [5 6 7 12 13 14];

    A_F_game = [A_F(:, colsPlayer_645), A_F(:, colsPlayer_611), ...
                A_F(:, colsPlayer_652), A_F(:, colsPlayer_671)];

    fprintf('[%s] A_F_game size: %d x %d\n', label, size(A_F_game,1), size(A_F_game,2));

    m_sh = 2*nPhys;
    assert(isequal(size(A_F_game), [m_sh, 14]), 'A_F_game has unexpected dimensions.');
    assert(isequal(size(b_shF), [m_sh, 1]), 'b_shF has unexpected dimensions.');
    assert(norm(A_F_game(nPhys+1:end,:) + A_F_game(1:nPhys,:), 'fro') < 1e-10, ...
        'Upper/lower voltage rows are inconsistent.');

    % Congestion diagnostics
    margin_upper = vmax - vstar_phys;
    margin_lower = vstar_phys - vmin;
    assert(norm(b_shF - [margin_upper; margin_lower], inf) < 1e-9, ...
        'b_shF does not match [vmax - vstar_phys; vstar_phys - vmin].');

    physicalFull = idx.nonSlackPhase(valid_ns).';
    busOf   = ceil(physicalFull/3);
    phaseOf = mod(physicalFull-1,3) + 1;
    playerBuses = [5 7 12 13];
    tol_small = 0.02;

    n_small_up = nnz(margin_upper < tol_small);
    n_small_lo = nnz(margin_lower < tol_small);
    n_neg_up   = nnz(margin_upper < 0);
    n_neg_lo   = nnz(margin_lower < 0);

    fprintf('[%s] load_scale=%.3f  vstar_phys range [%.4f,%.4f]  small-margin rows up/lo %d/%d  violated up/lo %d/%d\n', ...
        label, load_scale, min(vstar_phys), max(vstar_phys), n_small_up, n_small_lo, n_neg_up, n_neg_lo);

    out.idx = idx;
    out.S_v = S_v;
    out.S_theta = S_theta;
    out.m_v = m_v;
    out.valid_ns = valid_ns;
    out.nPhys = nPhys;
    out.vstar_phys = vstar_phys;
    out.vmin = vmin;
    out.vmax = vmax;
    out.A_grid = A_grid;
    out.b_grid = b_grid;
    out.A_F_game = A_F_game;
    out.b_shF = b_shF;
    out.m_sh = m_sh;
    out.margin_upper = margin_upper;
    out.margin_lower = margin_lower;
    out.busOf = busOf;
    out.phaseOf = phaseOf;
    out.player_ids = [645; 611; 652; 671];
    out.n_physical = [4; 2; 2; 6];
    out.block_start = [1; 5; 7; 9];
    out.block_end   = [4; 6; 8; 14];
    out.coordinate_convention = "physical variables are active/reactive deviations from baseline";
    out.constraint_convention = "A_F_game * Delta_s_F_game <= b_shF";
end


function export_shared_constraint(fit, filename, load_scale, relin_meta)
%EXPORT_SHARED_CONSTRAINT Save a build_shared_constraint() result in the
%same layout/variable names as example_ieee13.m's export, plus relin_meta
%(empty struct() for the flat model, iteration diagnostics for the
%congested-point model).

    A_F_game = fit.A_F_game;
    b_shF = fit.b_shF;
    player_ids = fit.player_ids;
    n_physical = fit.n_physical;
    block_start = fit.block_start;
    block_end = fit.block_end;
    coordinate_convention = fit.coordinate_convention;
    constraint_convention = fit.constraint_convention;

    stress_meta = struct( ...
        'load_scale', load_scale, ...
        'vmin',       fit.vmin(1), ...
        'vmax',       fit.vmax(1), ...
        'vstar_min',  min(fit.vstar_phys), ...
        'vstar_max',  max(fit.vstar_phys), ...
        'm_sh',       fit.m_sh, ...
        'relin',      relin_meta);

    export_dir = fullfile(pwd, "data");
    if ~exist(export_dir, "dir")
        mkdir(export_dir);
    end
    export_file = fullfile(export_dir, filename);

    save(export_file, ...
        "A_F_game", "b_shF", "player_ids", "n_physical", ...
        "block_start", "block_end", "coordinate_convention", ...
        "constraint_convention", "stress_meta", "-v7");

    fprintf('\nExported network data to:\n%s\n', export_file);
end
