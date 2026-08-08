function [v_lin, t_lin, Amat, Bmat] = bd_linear_solve(Y, n, p, q, v_slack, t_slack)
%BD_LINEAR_SOLVE One-shot Bolognani-Dorfler / 1ACPF manifold solve.
%
%   [v_lin, t_lin, Amat, Bmat] = BD_LINEAR_SOLVE(Y, n, p, q, v_slack, t_slack)
%
%   This is exactly the method of the original example_ieee13.m / the
%   Allerton 2015 paper: implicit linearization of the power-flow manifold
%   around the FLAT / no-load reference voltage profile (every bus-phase at
%   the ideal balanced triple aaa = [1; a; a^2]), giving a single linear
%   system whose solution approximates the true nonlinear operating point
%   in one solve -- no Newton iteration. Only the right-hand side (the
%   specified p,q and slack v,theta) changes from one call to the next; the
%   left-hand operator (Amat's top-left block) does not depend on loading.
%
%   This is "the manifold solve" (steps 1-2 of each stage in the DSO
%   pipeline): it produces x* = (v*, theta*, p*, q*), the approximate
%   operating point that stage 3 (bd_tangent_general.m) then re-linearizes
%   around exactly, and that is independently checked against a nonlinear
%   AC solve.
%
%   INPUTS
%   p, q             : 3n x 1 net-consumption specification at all
%                       non-slack bus-phases (values at slack rows unused).
%   v_slack, t_slack : 3 x 1 slack-bus voltage magnitude / angle.
%
%   OUTPUTS
%   v_lin, t_lin : 3n x 1 manifold-solution voltage magnitude / angle.
%   Amat, Bmat   : the assembled linear system, x = Amat \ Bmat, returned
%                  for diagnostics (residuals, conditioning) if needed.

    e0  = [1; zeros(n-1,1)];
    a   = exp(-1j*2*pi/3);
    aaa = [1; a; a^2];

    VTV = [kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n), zeros(3, 3*n)];
    VTT = [zeros(3, 3*n), kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n)];
    PQP = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1)), zeros(3*(n-1),3*n)];
    PQQ = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1))];

    NNN = Nmatrix(6*n);
    LLL = bracket(Y);
    RRR = bracket(kron(eye(n), diag(aaa)));   % flat/no-load reference

    Amat = [NNN*inv(RRR)*LLL*RRR, eye(6*n); VTV; VTT; PQP; PQQ];
    Bmat = [zeros(3*n,1); zeros(3*n,1); v_slack; t_slack; p(rw(2:n)); q(rw(2:n))];

    x = Amat \ Bmat;

    v_lin = x(1:3*n);
    t_lin = x(3*n+1:2*3*n);
end
