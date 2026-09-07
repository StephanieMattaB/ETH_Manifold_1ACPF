function res = min_effort_qp(A_sh, b_sh, A_loc, b_loc, opts)
%MIN_EFFORT_QP Minimum-effort corrective dispatch:
%
%       minimize_x  (1/2) * x'x
%       subject to  A_sh*x  <= b_sh
%                   A_loc*x <= b_loc
%
%   res = MIN_EFFORT_QP(A_sh, b_sh, A_loc, b_loc, opts)
%
%   No quadprog/optim package is available in this environment (no network
%   access to install one), so this is solved via ADMM -- a standard,
%   provably-convergent method for exactly this problem shape (strongly
%   convex quadratic objective, polyhedral constraint set):
%
%       min_x  (1/2)x'x + g(z)   s.t.  Cx - z = 0,   g(z) = indicator{z<=d}
%
%   x-update is a single small (nF x nF) linear solve; z-update is a
%   closed-form elementwise projection (min(.,d)); both are exact, so the
%   only approximation is the number of ADMM iterations run. The result is
%   NOT trusted blindly: after iterating, this function verifies the
%   solution against the constraints directly (constraint residual) rather
%   than relying on ADMM's own convergence diagnostics alone.
%
%   opts.rho (default 2), opts.max_iter (default 20000), opts.tol (1e-10).

    if nargin < 5; opts = struct(); end
    if ~isfield(opts,'rho');      opts.rho = 2; end
    if ~isfield(opts,'max_iter'); opts.max_iter = 20000; end
    if ~isfield(opts,'tol');      opts.tol = 1e-10; end

    C = [A_sh; A_loc];
    d = [b_sh; b_loc];
    [m, nF] = size(C);
    rho = opts.rho;

    % Cache the factorization of (I + rho*C'C), fixed across all iterations.
    M = eye(nF) + rho*(C'*C);
    R = chol(M);   % M = R'R, solve via two triangular solves each iteration

    z = zeros(m,1);
    u = zeros(m,1);
    x = zeros(nF,1);

    for k = 1:opts.max_iter
        rhs = rho*C'*(z - u);
        x = R \ (R' \ rhs);

        Cx = C*x;
        z_new = min(Cx + u, d);
        u = u + Cx - z_new;

        r_norm = norm(Cx - z_new, inf);           % primal residual
        s_norm = rho*norm(C'*(z_new - z), inf);    % dual residual
        z = z_new;

        if r_norm < opts.tol && s_norm < opts.tol
            break;
        end
    end

    res.x = x;
    res.iterations = k;
    res.converged = (r_norm < 1e-6 && s_norm < 1e-6);
    res.primal_residual = r_norm;
    res.dual_residual = s_norm;

    % Independent verification: does x actually satisfy the constraints,
    % regardless of what ADMM's own residuals claim?
    res.max_constraint_violation = max(C*x - d);
    res.feasible = res.max_constraint_violation <= 1e-6;
    if ~res.feasible
        warning('min_effort_qp:notFeasible', ...
            'ADMM result violates constraints by %.3e after %d iterations -- not trustworthy.', ...
            res.max_constraint_violation, k);
    end

    res.norm2 = norm(x);
end
