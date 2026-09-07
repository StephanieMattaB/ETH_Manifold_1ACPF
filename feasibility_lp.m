function res = feasibility_lp(A_sh, b_sh, A_loc, b_loc)
%FEASIBILITY_LP Pure feasibility test for F = {x : A_sh*x <= b_sh, A_loc*x <= b_loc}.
%
%   res = FEASIBILITY_LP(A_sh, b_sh, A_loc, b_loc)
%
%   Solves the standard slack-minimization feasibility LP via glpk:
%
%       min_{x,t}  t
%       s.t.       A_sh*x - t <= b_sh     (shared voltage-security rows get slack)
%                  A_loc*x    <= b_loc    (local device limits: HARD, no slack)
%                  t free
%
%   t* <= 0  =>  F is non-empty (a feasible x exists; t*<0 means strict
%                interior slack on the worst shared row).
%   t* >  0  =>  F is empty; t* is the minimum unavoidable violation of the
%                worst shared row, even in the best case.
%
%   Local (device) constraints are never relaxed -- they are physical
%   limits, not something a slack variable should be allowed to violate.

    nF = size(A_sh, 2);
    m_sh = size(A_sh, 1);
    m_loc = size(A_loc, 1);

    % variables z = [x (nF); t (1)]
    c = [zeros(nF,1); 1];

    A = [A_sh, -ones(m_sh,1); A_loc, zeros(m_loc,1)];
    b = [b_sh; b_loc];
    ctype = repmat('U', m_sh + m_loc, 1);   % all rows are <=

    lb = [-1e3*ones(nF,1); -1e3];
    ub = [ 1e3*ones(nF,1);  1e3];
    vartype = repmat('C', nF+1, 1);

    % NOTE: with only 3 output args, glpk's 3rd return is ERRNUM (0 = solver
    % invocation succeeded), NOT the solution status -- the actual GLPK
    % status code (5 = optimal) is in the 4th output, extra.status. Asking
    % for 3 outputs and treating that as a status code silently misreads
    % errnum==0 as "not optimal" (status codes never equal 0), so request
    % all 4 explicitly.
    [zopt, tmin, errnum, extra] = glpk(c, A, b, lb, ub, ctype, vartype, 1);

    res.errnum = errnum;
    res.status = extra.status;   % glpk: 5 = optimal
    res.solved = (errnum == 0) && (extra.status == 5);
    if ~res.solved
        error('feasibility_lp:solverFailed', ...
            'glpk did not return an optimal solution (errnum=%d, status=%d).', errnum, extra.status);
    end

    res.x = zopt(1:nF);
    res.t = tmin;
    res.feasible = (tmin <= 1e-9);

    res.slack_sh  = b_sh  - A_sh*res.x;
    res.slack_loc = b_loc - A_loc*res.x;
end
