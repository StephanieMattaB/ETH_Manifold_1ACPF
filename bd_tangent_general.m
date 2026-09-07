function tp = bd_tangent_general(Y, n, v_star, t_star)
%BD_TANGENT_GENERAL Local tangent-plane linearization of the power-flow
%manifold at an ARBITRARY operating point (v_star, t_star), i.e. the
%general (non-flat) Bolognani-Dorfler formulation (their eq. 4-5),
%specialized to the flat/no-load point ONLY as a sanity check, never as
%the working formula.
%
%   tp = BD_TANGENT_GENERAL(Y, n, v_star, t_star)
%
%   This is "step 3" of every stage in the DSO pipeline: after the manifold
%   solve (bd_linear_solve.m) produces x* = (v*, theta*, ...), this function
%   builds the exact tangent plane of the manifold AT that point:
%
%       A_x* = [ M(u*) * R(u*),  I ]
%       M(u*) = bracket(diag(conj(Y*u*))) + bracket(diag(u*))*N*bracket(Y)
%       R(u*) = Rmatrix(v*, theta*)
%
%   where p,q follow this repository's NET-CONSUMPTION sign convention
%   (positive = drawn from the network), matching ieee13.m / bd_linear_solve.
%
%   A fresh call to this function at a new (v_star, t_star) is REQUIRED at
%   every stage of the pipeline (S, S', S'') -- the linearization from a
%   previous stage must never be reused, since M(u*) and R(u*) both depend
%   on the point they are evaluated at.
%
%   OUTPUT (struct tp)
%   Atop     : 6n x 6n block  M(u*)*R(u*)  (the LHS operator on [v;theta]).
%   Amat     : full 12n x 12n assembled system [Atop, I; VTV;VTT;PQP;PQQ],
%              same bus-model rows (slack pinned, non-slack p,q selected)
%              as bd_linear_solve.m, for residual/conditioning checks.
%   Sstar    : 3n x 1 complex self-consistent injection u*.*conj(Y*u*) at
%              this point (the TRUE nonlinear injection the manifold would
%              need to be exactly on-manifold here) -- comparing -Sstar
%              against the specified net consumption is a residual
%              diagnostic of how far x* sits from the true manifold.
%   rcond_top: reciprocal condition number of Atop (conditioning diagnostic).

    a   = exp(-1j*2*pi/3);
    aaa = [1; a; a^2];
    e0  = [1; zeros(n-1,1)];

    u_star = v_star .* exp(1j*t_star);

    NNN = Nmatrix(6*n);
    Mmat = bracket(diag(conj(Y*u_star))) + bracket(diag(u_star))*NNN*bracket(Y);
    Rmat = Rmatrix(v_star, t_star);
    Atop = Mmat * Rmat;

    % One-time self-check: at the EXACT flat point, this must reduce to
    % bd_linear_solve.m's flat-only operator, to machine precision.
    if nargin == 4 && isequal(v_star, ones(3*n,1)) && isequal(t_star, kron(ones(n,1), angle(aaa)))
        RRR_flat = bracket(kron(eye(n), diag(aaa)));
        Atop_flat_ref = NNN*inv(RRR_flat)*bracket(Y)*RRR_flat;
        rel = norm(Atop - Atop_flat_ref, 'fro') / norm(Atop_flat_ref, 'fro');
        if rel > 1e-8
            error('bd_tangent_general:selfcheck', ...
                'General formula does not reduce to the flat-only operator at the flat point (rel diff %.3e).', rel);
        end
    end

    VTV = [kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n), zeros(3, 3*n)];
    VTT = [zeros(3, 3*n), kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n)];
    PQP = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1)), zeros(3*(n-1),3*n)];
    PQQ = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1))];

    tp.Atop  = Atop;
    tp.Amat  = [Atop, eye(6*n); VTV; VTT; PQP; PQQ];
    tp.Sstar = u_star .* conj(Y*u_star);
    % NOTE: Atop alone (no slack pinning) is structurally singular, exactly
    % like a bus admittance matrix without a reference bus -- that is
    % expected and not a conditioning problem. The meaningful conditioning
    % diagnostic is the FULL assembled (square, invertible) system.
    tp.rcond_top = rcond(tp.Amat);
    tp.v_star = v_star;
    tp.t_star = t_star;
end
