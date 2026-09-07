function results = validate_io_map(map, sh, physics, opts)
% VALIDATE_IO_MAP  Run the four map-construction tests from section 12 of the
%                  bridge doc, BEFORE any GNE solve. Cheap bug-catchers.
%
%   results = validate_io_map(map, sh, physics, opts)
%
%   map      : struct from build_io_map
%   sh       : struct from build_shared_constraint (may be [] to skip 3 & 4)
%   physics  : struct with fields Ybus, slack, Vslack for the finite-difference
%              test (Test 2). If [], Test 2 is skipped.
%   opts     : struct, optional. opts.eps (FD step, default 1e-6),
%              opts.tol_fd (default 1e-4), opts.verbose (default true).
%
%   Returns a struct of pass/fail flags and residuals; prints a summary.

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'eps'),     opts.eps = 1e-6;   end
    if ~isfield(opts,'tol_fd'),  opts.tol_fd = 1e-4; end
    if ~isfield(opts,'verbose'), opts.verbose = true; end

    S = map.S; m = map.m;
    results = struct();

    % ---- Test 1: baseline consistency  ||y* - (S u* + m)||inf ~ 0 -----------
    r1 = norm(map.y_star - (S*map.u_star + m), inf);
    results.test1_resid = r1;
    results.test1_pass  = r1 < 1e-9;

    % ---- Test 2: finite-difference sensitivity vs columns of S -------------
    results.test2_pass = NaN; results.test2_maxerr = NaN;
    if ~isempty(physics)
        Ybus = physics.Ybus; slack = physics.slack; Vslack = physics.Vslack;
        n = map.u_idx.n; nonslack = map.nonslack;
        % reconstruct baseline injections (absolute) from u_star
        n_ns = numel(nonslack);
        Pbase = zeros(n,1); Qbase = zeros(n,1);
        Pbase(nonslack) = map.u_star(1:n_ns);
        Qbase(nonslack) = map.u_star(n_ns+1:end);
        % baseline AC solve (warm start from lin. point)
        [Vm0,Va0] = nrpf(Ybus,Pbase,Qbase,slack,Vslack);
        y0 = Vm0(map.mon_buses);

        Scols_fd = zeros(size(S));
        ep = opts.eps;
        for k = 1:size(S,2)
            Pk = Pbase; Qk = Qbase;
            if k <= n_ns, Pk(nonslack(k)) = Pk(nonslack(k)) + ep;
            else,         Qk(nonslack(k-n_ns)) = Qk(nonslack(k-n_ns)) + ep; end
            [Vmk,~] = nrpf(Ybus,Pk,Qk,slack,Vslack,Vm0,Va0);
            yk = Vmk(map.mon_buses);
            Scols_fd(:,k) = (yk - y0)/ep;
        end
        err2 = max(abs(Scols_fd(:) - S(:)));
        results.test2_maxerr = err2;
        results.test2_pass = err2 < opts.tol_fd;
        results.S_fd = Scols_fd;
    end

    % ---- Test 3: protected-load absorption ---------------------------------
    % A_F * u_hat_F <= b_sh'   <=>   A_grid*[u_hat_F;u_hat_P] <= b_grid
    results.test3_pass = NaN; results.test3_resid = NaN;
    if ~isempty(sh)
        u_hat_F = map.u_star(sh.FL_cols);
        lhs_partition = sh.A_F * u_hat_F;              % uses b_sh_prime
        % full-grid check with the same baseline
        lhs_grid = sh.A_sh_grid * map.u_star;
        slack3   = lhs_partition - sh.b_sh_prime;
        slack3g  = lhs_grid      - sh.b_sh_grid;
        r3 = norm(slack3 - slack3g, inf);
        results.test3_resid = r3;
        results.test3_pass  = r3 < 1e-9;
    end

    % ---- Test 4: zero-gamma physics ----------------------------------------
    % third column of every per-player block must be identically zero
    results.test4_pass = NaN; results.test4_maxabs = NaN;
    if ~isempty(sh)
        mx = 0;
        for j = 1:numel(sh.A_sh_blocks)
            mx = max(mx, max(abs(sh.A_sh_blocks{j}(:,3))));
        end
        results.test4_maxabs = mx;
        results.test4_pass = (mx == 0);
    end

    if opts.verbose
        fprintf('\n===== validate_io_map =====\n');
        pf = @(b) merge(isnan(b),'  skip', merge(logical(b),'  PASS','  FAIL'));
        fprintf('Test 1 baseline consistency : %s  (resid %.2e)\n', pf(results.test1_pass), r1);
        fprintf('Test 2 finite-difference    : %s  (max err %.2e)\n', pf(results.test2_pass), results.test2_maxerr);
        fprintf('Test 3 protected absorption : %s  (resid %.2e)\n', pf(results.test3_pass), results.test3_resid);
        fprintf('Test 4 zero-gamma physics   : %s  (max |col3| %.2e)\n', pf(results.test4_pass), results.test4_maxabs);
        fprintf('===========================\n');
    end
end

function s = merge(cond, a, b)
    if cond, s = a; else, s = b; end
end
