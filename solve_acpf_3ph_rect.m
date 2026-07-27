function result = solve_acpf_3ph_rect(Y, slack_idx, Vslack, Sspec, V0, opts)
%SOLVE_ACPF_3PH_RECT Nonlinear three-phase AC power flow in rectangular form.
%
% Solves, at every non-slack phase terminal,
%
%       S_calc(V) = V .* conj(Y*V) = S_spec
%
% where S_spec follows the network-injection convention:
%   positive P/Q = injection into the network,
%   negative P/Q = consumption from the network.
%
% INPUTS
%   Y          N-by-N complex nodal admittance matrix [p.u.]
%   slack_idx  indices of fixed-voltage phase terminals
%   Vslack     complex slack voltages corresponding to slack_idx [p.u.]
%   Sspec      N-by-1 specified complex network injections [p.u.]
%   V0         optional N-by-1 initial complex voltage vector [p.u.]
%   opts       optional structure:
%                .tol               default 1e-10
%                .max_iter          default 50
%                .max_line_search   default 15
%                .verbose           default true
%
% OUTPUT
%   result.V                 complex terminal voltages
%   result.v                 voltage magnitudes
%   result.theta             voltage angles [rad]
%   result.S_calc            calculated complex injections
%   result.S_spec            specified complex injections
%   result.slack_power       calculated slack complex powers
%   result.mismatch          complex non-slack mismatch
%   result.max_mismatch      infinity norm of rectangular mismatch
%   result.converged         logical convergence flag
%   result.iterations        Newton iterations
%   result.mismatch_history  residual history
%   result.step_history      accepted step-size history
%
% The solver assumes fixed-PQ injections at all non-slack terminals. This is
% the appropriate first nonlinear validation for the current 1ACPF/GNE model,
% in which p and q are explicit fixed/controlled quantities.

    if nargin < 4
        error('At least Y, slack_idx, Vslack, and Sspec are required.');
    end
    if nargin < 5
        V0 = [];
    end
    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'tol'),             opts.tol = 1e-10; end
    if ~isfield(opts, 'max_iter'),        opts.max_iter = 50; end
    if ~isfield(opts, 'max_line_search'), opts.max_line_search = 15; end
    if ~isfield(opts, 'verbose'),         opts.verbose = true; end

    [N, M] = size(Y);
    if N ~= M
        error('Y must be square.');
    end

    slack_idx = slack_idx(:);
    Vslack = Vslack(:);
    Sspec = Sspec(:);

    if numel(Vslack) ~= numel(slack_idx)
        error('Vslack and slack_idx must have the same length.');
    end
    if numel(Sspec) ~= N
        error('Sspec must have one entry per terminal.');
    end
    if any(~isfinite(real(Y(:)))) || any(~isfinite(imag(Y(:))))
        error('Y contains NaN or Inf.');
    end

    all_idx = (1:N).';
    ns_idx = setdiff(all_idx, slack_idx, 'stable');
    n_ns = numel(ns_idx);

    if isempty(V0)
        n_slack = numel(slack_idx);
        if mod(N, n_slack) ~= 0
            error(['Automatic initialization requires the number of terminals ', ...
                   'to be a multiple of the number of slack phases.']);
        end
        V = repmat(Vslack, N/n_slack, 1);
    else
        V = V0(:);
        if numel(V) ~= N
            error('V0 must have one entry per terminal.');
        end
    end

    V(slack_idx) = Vslack;

    G = real(Y);
    B = imag(Y);

    mismatch_history = NaN(opts.max_iter + 1, 1);
    step_history = NaN(opts.max_iter, 1);
    converged = false;
    iterations = 0;

    for k = 0:opts.max_iter
        [F, Scalc] = mismatch_vector(Y, V, Sspec, ns_idx);
        mismatch_norm = norm(F, inf);
        mismatch_history(k + 1) = mismatch_norm;

        if opts.verbose
            fprintf('Newton iteration %2d: ||F||_inf = %.3e\n', ...
                k, mismatch_norm);
        end

        if mismatch_norm <= opts.tol
            converged = true;
            iterations = k;
            break;
        end

        if k == opts.max_iter
            iterations = k;
            break;
        end

        e = real(V);
        f = imag(V);

        Ire = G*e - B*f;
        Iim = B*e + G*f;

        % Analytic Jacobian of [P;Q] with respect to [e;f].
        JPe = diag(Ire) + diag(e)*G + diag(f)*B;
        JPf = diag(Iim) - diag(e)*B + diag(f)*G;
        JQe = diag(f)*G - diag(Iim) - diag(e)*B;
        JQf = diag(Ire) - diag(f)*B - diag(e)*G;

        J = [ ...
            JPe(ns_idx, ns_idx), JPf(ns_idx, ns_idx); ...
            JQe(ns_idx, ns_idx), JQf(ns_idx, ns_idx) ...
        ];

        if any(~isfinite(J(:)))
            error('The Newton Jacobian contains NaN or Inf.');
        end

        if rcond(full(J)) < 1e-14
            warning('The Newton Jacobian is nearly singular at iteration %d.', k);
        end

        dz = -J \ F;
        dV = dz(1:n_ns) + 1i*dz(n_ns + 1:end);

        % Backtracking line search. It improves robustness for stressed cases.
        alpha = 1.0;
        accepted = false;
        best_V = V;
        best_norm = mismatch_norm;

        for ls = 1:opts.max_line_search
            Vtrial = V;
            Vtrial(ns_idx) = V(ns_idx) + alpha*dV;
            Vtrial(slack_idx) = Vslack;

            Ftrial = mismatch_vector(Y, Vtrial, Sspec, ns_idx);
            trial_norm = norm(Ftrial, inf);

            if trial_norm < best_norm
                best_norm = trial_norm;
                best_V = Vtrial;
            end

            % Armijo-type sufficient decrease.
            if trial_norm <= (1 - 1e-4*alpha)*mismatch_norm
                V = Vtrial;
                accepted = true;
                break;
            end

            alpha = alpha/2;
        end

        if ~accepted
            % Use the best tested damped step rather than taking a blind full step.
            if best_norm < mismatch_norm
                V = best_V;
            else
                warning(['Newton step failed to reduce the residual at iteration %d. ', ...
                         'Stopping without declaring convergence.'], k);
                iterations = k;
                break;
            end
        end

        step_history(k + 1) = alpha;
        iterations = k + 1;
    end

    [Ffinal, Scalc] = mismatch_vector(Y, V, Sspec, ns_idx);
    complex_mismatch = Scalc(ns_idx) - Sspec(ns_idx);

    result = struct();
    result.V = V;
    result.v = abs(V);
    result.theta = angle(V);
    result.S_calc = Scalc;
    result.S_spec = Sspec;
    result.slack_power = Scalc(slack_idx);
    result.mismatch = complex_mismatch;
    result.max_mismatch = norm(Ffinal, inf);
    result.converged = converged;
    result.iterations = iterations;
    result.slack_idx = slack_idx;
    result.non_slack_idx = ns_idx;
    result.mismatch_history = mismatch_history(1:iterations + 1);

    if iterations > 0
        result.step_history = step_history(1:iterations);
    else
        result.step_history = [];
    end
end


function [F, Scalc] = mismatch_vector(Y, V, Sspec, ns_idx)
%MISMATCH_VECTOR Rectangular P/Q mismatch at non-slack terminals.

    Scalc = V .* conj(Y*V);
    mismatch = Scalc(ns_idx) - Sspec(ns_idx);

    F = [real(mismatch); imag(mismatch)];
end
