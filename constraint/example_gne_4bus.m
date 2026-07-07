% EXAMPLE_GNE_4BUS  End-to-end demo of the ACPF-manifold -> shared-constraint
% bridge on a small MESHED 4-bus feeder, with NO MatPower dependency.
%
% Pipeline (section 14 of the bridge doc):
%   1-2  build Ybus, get baseline AC operating point (nrpf)
%   3-9  build tangent map -> S, m
%   10   grid-level shared constraint
%   11-12 validation Tests 1-2 (baseline, finite-difference)
%   13   FL/protected partition -> b_sh'
%   14   per-player Pi blocks -> A_sh blocks ; validation Tests 3-4

clear; clc;

% ---- 1. grid: 4-bus mesh -------------------------------------------------
%  bus 1 = slack, buses 2,3 = FL (flexible), bus 4 = protected load
n = 4; slack = 1; Vslack = 1.0;
lines = [ ...  from  to     R       X
            1   2   0.02   0.06 ;
            2   3   0.03   0.09 ;
            3   4   0.025  0.07 ;
            1   4   0.02   0.05 ;
            2   4   0.04   0.08 ];
Ybus = zeros(n);
for l = 1:size(lines,1)
    i = lines(l,1); j = lines(l,2);
    y = 1/(lines(l,3) + 1j*lines(l,4));
    Ybus(i,i) = Ybus(i,i) + y;  Ybus(j,j) = Ybus(j,j) + y;
    Ybus(i,j) = Ybus(i,j) - y;  Ybus(j,i) = Ybus(j,i) - y;
end

% specified injections (loads negative). Slack entry ignored.
Pspec = [ 0 ; -0.30 ; -0.20 ; -0.25 ];
Qspec = [ 0 ; -0.10 ; -0.05 ; -0.08 ];

% ---- 2. baseline AC operating point -> linearization point ----------------
[V0, A0, conv] = nrpf(Ybus, Pspec, Qspec, slack, Vslack);
assert(conv, 'baseline AC power flow did not converge');
fprintf('Baseline AC solution:\n');
for b = 1:n, fprintf('  bus %d: |V| = %.5f, angle = %+.4f deg\n', b, V0(b), A0(b)*180/pi); end

% ---- 3-9. build the input-output map --------------------------------------
mon_buses = [2 3];                 % monitor voltage magnitude at buses 2 and 3
map = build_io_map(Ybus, V0, A0, slack, mon_buses);
fprintf('\nSensitivity S (rows = |V| at buses %s ; cols = [P2 P3 P4 Q2 Q3 Q4]):\n', mat2str(mon_buses));
disp(map.S);
fprintf('offset m = %s ,  rcond(A_U) = %.2e\n', mat2str(map.m',4), map.rcond_AU);

% ---- 10-13. shared constraint, partition, per-player blocks ---------------
y_lo = [0.95; 0.95];  y_hi = [1.05; 1.05];
FL_buses = [2 3];                  % strategic flexible players; bus 4 protected
sh = build_shared_constraint(map, y_lo, y_hi, FL_buses);
fprintf('\nb_sh'' (protected load absorbed) = %s\n', mat2str(sh.b_sh_prime',4));
fprintf('per-player block A_sh{1} (bus 2), 3rd col should be 0:\n'); disp(sh.A_sh_blocks{1});

% ---- 11-12 & 3-4. validation ---------------------------------------------
physics = struct('Ybus',Ybus,'slack',slack,'Vslack',Vslack);
res = validate_io_map(map, sh, physics);

fprintf('\nAll four tests passed: %d\n', ...
    res.test1_pass && res.test2_pass && res.test3_pass && res.test4_pass);
