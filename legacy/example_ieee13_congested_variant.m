%	This code generates the simulations included in
%
%	S. Bolognani, F. Dörfler (2015)
%	"Fast power system analysis via implicit linearization of the power flow manifold."
%	In Proc. 53rd Annual Allerton Conference on Communication, Control, and Computing.
%	Preprint available at http://control.ee.ethz.ch/~bsaverio/papers/BolognaniDorfler_Allerton2015.pdf
%
%	This source code is distributed in the hope that it will be useful, but without any warranty.
%
%	MatLab OR GNU Octave, version 3.8.1 available at http://www.gnu.org/software/octave/
%	MATPOWER 5.1 available at http://www.pserc.cornell.edu/matpower/
%
%   Commented version by Stephanie Matta Spring 2026 KTH
%   File solves and compare, linked voltage profiles - network dependent
%   Helpers dependent (Nmatrix, Rmatrix, rw and bracket.m)
%   Mx=b


%clear all 
close all

Vbase = 4160/sqrt(3); %phase-to-line base voltage
Sbase = 5e6;          %VA
Zbase = Vbase^2/Sbase;%Ohms

[Y,v,t,p,q,n] = ieee13();

% initial voltage profiles 3x13 hard-coded
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

v_testfeeder = reshape(v_testfeeder.',3*n,1); %to col vector node-phase 3nx1

%begins
t_testfeeder = [...
    0.00    -120.00 120.00  ;...
    -2.49   -121.72 117.83  ;...
    -2.56   -121.77 117.82  ;...
    -3.23   -122.22 117.34  ;...
    NaN     -121.90 117.86  ;...
    NaN     -121.98 117.90  ;...
    -5.30   -122.34 116.02  ;...
    -5.31   -122.34 116.02  ;...
    -5.56   -122.52 116.03  ;...
    -5.30   -122.34 116.02  ;...
    -5.32   NaN     115.92  ;...
    NaN     NaN     115.78  ;...
    -5.25   NaN     NaN     ];

t_testfeeder = reshape(t_testfeeder.',3*n,1)/180*pi; %col vector (&degrees-rad)

%%

% Linearized model

e0 = [1;zeros(n-1,1)]; %selecting slack bus (bus 1)
a = exp(-1j*2*pi/3);   %phasor 120d rotation
aaa = [1; a; a^2];     %phase sequence vector 

%constraint matrix enforcing vector magnitudes
VTV = [kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n), zeros(3, 3*n)];
%constraint matrix enforcing vector angles
VTT = [zeros(3, 3*n), kron(e0',eye(3)), zeros(3, 3*n), zeros(3, 3*n)];
%block selecting P unkowns for non-slack-bus
PQP = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1)), zeros(3*(n-1),3*n)];
%block selecting Q unkowns for non-slack-bus
PQQ = [zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3*n), zeros(3*(n-1),3), eye(3*(n-1))];
%blkdiag phase-rotation + helper
UUU = bracket(kron(eye(n),diag(aaa)));
%indexing selection matrix N for 6n-size linear sys
NNN = Nmatrix(6*n);
%bracked admitance matrix
LLL = bracket(Y);
%rotation R for phase angles + helper
PPP = Rmatrix(ones(3*n,1), kron(ones(n,1),angle(aaa)));


%% ============================================================
%  Congestion via operating-point stress (real limits kept)
%  ------------------------------------------------------------
%  Scale the non-slack demand so several operating-point voltages
%  sag toward the *real* limits vmin = 0.95 / vmax = 1.05.
%
%  NOTE on this model (no-load linearization):
%    * Amat is fixed => the sensitivity S_v (and therefore A_F_game)
%      is operating-point INDEPENDENT and does NOT change with load.
%      This is expected; a load-dependent S_v would need a full
%      nonlinear AC solve + Jacobian at the stressed point.
%    * The no-load offset m_v = vstar_ns - S_v*sstar_ns is invariant
%      to loading (verified below).
%    * What moves toward the limits is the operating-point voltage
%      vstar_ns, and the EXPORTED margins are exactly
%         b_shF = [vmax - vstar_ns ; vstar_ns - vmin].
%      So load scaling is the correct physical knob: it tightens b_shF.
%% ============================================================

load_scale = 1.8;   % global demand multiplier (sweep 1.5 .. 2.5) added

% Optional extra stress at electrically-weak / player-adjacent buses.
% Internal bus numbering; FaaS players sit at buses 5, 7, 12, 13.
% Add entries to drive specific feeder ends closer to a limit.
bus_extra = containers.Map('KeyType','double','ValueType','double');
% bus_extra(9)  = 1.3;   % heavy load deep in the feeder
% bus_extra(13) = 1.2;   % remote single-phase tap (player 652)

% handling correcly the capacitor scaling when stressing the network
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

% equivalent, when linearizing around the no load solution
RRR = bracket(kron(eye(n),diag(aaa))); %same as UUU
Amat = [NNN*inv(RRR)*LLL*RRR eye(6*n); VTV; VTT; PQP; PQQ]; %A Assembles linear sys
Bmat = [zeros(3*n,1); zeros(3*n,1); v(rw(1)); t(rw(1)); p(rw(2:n)); q(rw(2:n))]; %RHS b vector

x = Amat\Bmat; %solving Ax=B

v_linearized = x(1:3*n); %extracts first 3 linear magnitude and angles below
t_linearized = x(3*n+1:2*3*n); 


%%
%masking NaN for comparison
v_linearized(isnan(v_testfeeder))=NaN;
v_linearized_nl(isnan(v_testfeeder))=NaN;
t_linearized(isnan(t_testfeeder))=NaN;

% print

%loop over busese
for bus = 1:n

	fprintf('--- bus %i\t|  ', bus);
	vbus = v_linearized((bus-1)*3+1:(bus-1)*3+3); %magnitude
	tbus = t_linearized((bus-1)*3+1:(bus-1)*3+3); %angle
	fprintf(1, '%f at %+f  |  ', [vbus(1) tbus(1)/pi*180]);
	fprintf(1, '%f at %+f  |  ', [vbus(2) tbus(2)/pi*180]);
	fprintf(1, '%f at %f\n', [vbus(3) tbus(3)/pi*180]);

end

% fprintf('A (size %dx%d):\n', size(A,1), size(A,2));
% disp(full(A));
% 
% fprintf('x* (size %dx1):\n', numel(xstar));
% disp(xstar);
% 
% fprintf('b (size %dx1):\n', numel(b));
% disp(b);

% % show equation form and residual
% fprintf('Equation: A * x* = b\n');
% res = A*xstar - b;
% fprintf('Residual norm ||A*x* - b||_2 = %g\n', norm(res));


%adding:

%% ============================================================
%  Voltage sensitivity constraint from manifold
%  Output:
%      v_ns = S_v*s_ns + m_v
%      A_grid*s_ns <= b_grid
%% ============================================================

%% objects

A = Amat;
b = Bmat;
xstar = x;

%% dimensions and indexing

nBus   = 13;
nPhase = 3;
nPhi   = nBus*nPhase;       % 39 bus-phase slots

assert(isequal(size(A), [4*nPhi, 4*nPhi]), ...
    'A must be 156 x 156.');

assert(isequal(size(b), [4*nPhi, 1]), ...
    'b must be 156 x 1.');

assert(isequal(size(xstar), [4*nPhi, 1]), ...
    'xstar must be 156 x 1.');

% State blocks:
% x = [v; theta; p; q]

idx.v     = 1:nPhi;
idx.theta = nPhi + (1:nPhi);
idx.p     = 2*nPhi + (1:nPhi);
idx.q     = 3*nPhi + (1:nPhi);

% Bus 1 is the three-phase slack bus.
idx.slackPhase    = 1:3;
idx.nonSlackPhase = 4:nPhi;

idx.slack_v     = idx.v(idx.slackPhase);
idx.slack_theta = idx.theta(idx.slackPhase);

idx.nonSlack_v     = idx.v(idx.nonSlackPhase);
idx.nonSlack_theta = idx.theta(idx.nonSlackPhase);
idx.nonSlack_p     = idx.p(idx.nonSlackPhase);
idx.nonSlack_q     = idx.q(idx.nonSlackPhase);

nNonSlackPhi = numel(idx.nonSlackPhase);   % 36

%% Network manifold (ROWS)
% 1:39  -> active-polist wer equations
% 40:78 -> reactive-power equations

A_net = A(1:2*nPhi, :);
b_net = b(1:2*nPhi);

rowsP = 1:nPhi;
rowsQ = nPhi + (1:nPhi);

rowsNonSlack = [ ...
    rowsP(idx.nonSlackPhase), ...
    rowsQ(idx.nonSlackPhase) ...
];

A_red = A_net(rowsNonSlack, :);
b_red = b_net(rowsNonSlack);

%% Partition reduced variables

% z_ns = [v_ns; theta_ns]
colsZ = [ ...
    idx.nonSlack_v, ...
    idx.nonSlack_theta ...
];

% s_ns = [p_ns; q_ns]
colsS = [ ...
    idx.nonSlack_p, ...
    idx.nonSlack_q ...
];

% Fixed slack reference z_0 = [v_0; theta_0]
colsZ0 = [ ...
    idx.slack_v, ...
    idx.slack_theta ...
];

A_z = A_red(:, colsZ);
A_s = A_red(:, colsS);
A_0 = A_red(:, colsZ0);

%cheking
assert(isequal(size(A_z), [2*nNonSlackPhi, 2*nNonSlackPhi]), ...
    'Unexpected dimensions for A_z.');

assert(rank(A_z) == size(A_z,2), ...
    'A_z is rank deficient.');

if rcond(A_z) < 1e-8
    warning('A_z is poorly conditioned: rcond = %.3e.', rcond(A_z));
end

%% reference operating point

zstar_ns = [ ...
    xstar(idx.nonSlack_v); ...
    xstar(idx.nonSlack_theta) ...
];

sstar_ns = [ ...
    xstar(idx.nonSlack_p); ...
    xstar(idx.nonSlack_q) ...
];

zstar_0 = [ ...
    xstar(idx.slack_v); ...
    xstar(idx.slack_theta) ...
];

vstar_ns = xstar(idx.nonSlack_v);

%% sensitivity map

% Deviation map:
% A_z*Delta z_ns + A_s*Delta s_ns = 0
% Delta z_ns = S_zs*Delta s_ns

S_zs = -(A_z \ A_s);

S_v     = S_zs(1:nNonSlackPhi, :);
S_theta = S_zs(nNonSlackPhi+1:end, :);

% Active- and reactive-power sensitivity blocks

S_vP = S_v(:, 1:nNonSlackPhi);
S_vQ = S_v(:, nNonSlackPhi+1:end);

%% (abs) affine voltage map
% v_ns = S_v*s_ns + m_v

m_v = vstar_ns - S_v*sstar_ns;

%% adding
valid_ns = ~isnan(v_testfeeder(idx.nonSlackPhase));
nPhys = nnz(valid_ns);   % expected: 29

S_v_phys    = S_v(valid_ns,:);
m_v_phys    = m_v(valid_ns);
vstar_phys  = vstar_ns(valid_ns);

% Self-consistency: m_v is the no-load offset and MUST be invariant to
% load_scale (it should match the slack magnitude propagated through the
% network, ~1.05 range). If this drifts when you change load_scale, the
% partition is wrong.
fprintf('\nm_v range (no-load offset): [%.4f, %.4f] p.u.\n', ...
    min(m_v), max(m_v));

%% voltage constraint (set  bandwidth)

vmin = 0.95*ones(nPhys,1); % real regulatory limits (kept)
vmax = 1.05*ones(nPhys,1);

A_grid = [
     S_v_phys;
    -S_v_phys
];

b_grid = [
     vmax - m_v_phys;
    -vmin + m_v_phys
];


%% validation

fullResidual = norm(A*xstar - b, 2);

reducedResidual = norm( ...
    A_z*zstar_ns + ...
    A_s*sstar_ns + ...
    A_0*zstar_0 - b_red, 2);

mapResidual = norm( ...
    vstar_ns - (S_v*sstar_ns + m_v), inf);

constraintResidual = A_grid*sstar_ns - b_grid;

fprintf('\n--- Full sensitivity constraint ---\n');
fprintf('Full manifold residual:    %.3e\n', fullResidual);
fprintf('Reduced manifold residual: %.3e\n', reducedResidual);
fprintf('Affine-map residual:       %.3e\n', mapResidual);
fprintf('rcond(A_z):                %.3e\n', rcond(A_z));

fprintf('\nS_v size:     %d x %d\n', size(S_v,1), size(S_v,2));
fprintf('A_grid size:  %d x %d\n', size(A_grid,1), size(A_grid,2));
fprintf('b_grid size:  %d x %d\n', size(b_grid,1), size(b_grid,2));

fprintf('\nMaximum voltage: %.6f p.u.\n', max(vstar_phys));
fprintf('Minimum voltage: %.6f p.u.\n', min(vstar_phys));
fprintf('Maximum constraint violation: %.6e p.u.\n', ...
    max(constraintResidual));

%% defining relevant constraint for fleixble community only
% Flexible agents:
%   645 -> internal bus 5, phases b,c
%   611 -> internal bus 12, phase c
%   652 -> internal bus 13, phase a
%   671 -> internal bus 7, phases a,b,c
%
% Starting constraint:
%   A_grid * s_ns <= b_grid
% with
%   s_ns = [p_ns; q_ns]
%
% Final deviation constraint:
%   A_F * u_F <= b_shF
%
% where
%   u_F = Delta s_F

nNonSlackPhi = 36;

%% Phase and bus-phase indexing

phase.a = 1;
phase.b = 2;
phase.c = 3;

% Full phase-slot index in the 39-slot representation
busPhaseIdx = @(bus, ph) 3*(bus - 1) + ph;

% Convert full phase-slot index to the non-slack 36-slot representation
toNonSlackIdx = @(fullIdx) fullIdx - 3;

%% Flexible bus-phase locations

idxF_645 = toNonSlackIdx([ ...
    busPhaseIdx(5, phase.b), ...
    busPhaseIdx(5, phase.c) ...
]);

idxF_611 = toNonSlackIdx( ...
    busPhaseIdx(12, phase.c) ...
);

idxF_652 = toNonSlackIdx( ...
    busPhaseIdx(13, phase.a) ...
);

idxF_671 = toNonSlackIdx([ ...
    busPhaseIdx(7, phase.a), ...
    busPhaseIdx(7, phase.b), ...
    busPhaseIdx(7, phase.c) ...
]);

idxF_phase = [ ...
    idxF_645, ...
    idxF_611, ...
    idxF_652, ...
    idxF_671 ...
];

assert(numel(unique(idxF_phase)) == numel(idxF_phase), ...
    'Duplicate flexible bus-phase indices detected.');

%% Flexible and protected columns

% s_ns = [p_ns; q_ns]

idxF_p = idxF_phase;
idxF_q = nNonSlackPhi + idxF_phase;

cols_F = [idxF_p, idxF_q];

cols_all = 1:(2*nNonSlackPhi);
cols_P = setdiff(cols_all, cols_F, 'stable');

%% Baseline injection partition

sbar_ns = sstar_ns;

sbar_F = sbar_ns(cols_F);
sbar_P = sbar_ns(cols_P);

%% Grid-constraint partition

A_F = A_grid(:, cols_F);
A_P = A_grid(:, cols_P);

%% Flexible-community residual margin

% Protected injections remain fixed at baseline:
%
%   s_P = sbar_P
%
% Flexible injections are written as:
%
%   s_F = sbar_F + u_F

b_shF = b_grid ...
      - A_F*sbar_F ...
      - A_P*sbar_P;

% Equivalent compact expression
b_shF_check = b_grid - A_grid*sbar_ns;

assert(norm(b_shF - b_shF_check, inf) < 1e-10, ...
    'Flexible/protected partition is inconsistent.');

%% Final flexible-community constraint
%
%   A_F * u_F <= b_shF

fprintf('\nFlexible-community shared constraint\n');
fprintf('Flexible phase locations: %d\n', numel(idxF_phase));
fprintf('Flexible variables:       %d\n', numel(cols_F));
fprintf('Protected variables:      %d\n', numel(cols_P));
fprintf('A_F size:                 %d x %d\n', size(A_F,1), size(A_F,2));
fprintf('A_P size:                 %d x %d\n', size(A_P,1), size(A_P,2));
fprintf('b_shF size:               %d x %d\n', size(b_shF,1), size(b_shF,2));

%% ordering decisions per bus-phase <----
%% Player-wise shared-constraint blocks
%
% Physical flexible-variable ordering in A_F:
%
% [p645_b, p645_c, p611_c, p652_a, p671_a, p671_b, p671_c, ...
%  q645_b, q645_c, q611_c, q652_a, q671_a, q671_b, q671_c]

colsPlayer_645 = [1 2 8 9];
colsPlayer_611 = [3 10];
colsPlayer_652 = [4 11];
colsPlayer_671 = [5 6 7 12 13 14];

A_F_645 = A_F(:, colsPlayer_645);
A_F_611 = A_F(:, colsPlayer_611);
A_F_652 = A_F(:, colsPlayer_652);
A_F_671 = A_F(:, colsPlayer_671);

A_F_game = [ ...
    A_F_645, ...
    A_F_611, ...
    A_F_652, ...
    A_F_671 ...
];

%% Final shared constraint
%   A_F_game * z <= b_shF
% with
%   z = col(z_645, z_611, z_652, z_671)

fprintf('\nPlayer-wise shared constraint\n');
fprintf('A_F_645 size: %d x %d\n', size(A_F_645,1), size(A_F_645,2));
fprintf('A_F_611 size: %d x %d\n', size(A_F_611,1), size(A_F_611,2));
fprintf('A_F_652 size: %d x %d\n', size(A_F_652,1), size(A_F_652,2));
fprintf('A_F_671 size: %d x %d\n', size(A_F_671,1), size(A_F_671,2));
fprintf('A_F_game size:     %d x %d\n', size(A_F_game,1), size(A_F_game,2));


%% export 

%% ============================================================
% MATLAB -> Julia export: physical shared network constraint
%
% A_F_game * Delta_s_F_game <= b_shF
%% ============================================================

% Expected player-grouped physical column ordering:
%
% Player 645, columns 1:4
%   [dp_645_b, dp_645_c, dq_645_b, dq_645_c]
%
% Player 611, columns 5:6
%   [dp_611_c, dq_611_c]
%
% Player 652, columns 7:8
%   [dp_652_a, dq_652_a]
%
% Player 671, columns 9:14
%   [dp_671_a, dp_671_b, dp_671_c, ...
%    dq_671_a, dq_671_b, dq_671_c]

m_sh = 2*nPhys;   % 58 = physical voltage rows (upper + lower)

assert(isequal(size(A_F_game), [m_sh, 14]), ...
    'A_F_game has unexpected dimensions.');

assert(isequal(size(b_shF), [m_sh, 1]), ...
    'b_shF has unexpected dimensions.');

% Verify upper/lower voltage-row structure
assert(norm( ...
    A_F_game(nPhys+1:end,:) + A_F_game(1:nPhys,:), ...
    'fro') < 1e-10, ...
    'Upper/lower voltage rows are inconsistent.');

% Player metadata
player_ids = [645; 611; 652; 671];

n_physical = [4; 2; 2; 6];

block_start = [1; 5; 7; 9];
block_end   = [4; 6; 8; 14];

% Coordinate convention
coordinate_convention = ...
    "physical variables are active/reactive deviations from baseline";

constraint_convention = ...
    "A_F_game * Delta_s_F_game <= b_shF";


%% ============================================================
%  Congestion diagnostics on the exported shared constraint
%  ------------------------------------------------------------
%  b_shF encodes the operating-point margins on the physical rows:
%     rows        1:nPhys   ->  vmax - vstar_phys  (headroom to upper limit)
%     rows nPhys+1:2*nPhys  ->  vstar_phys - vmin  (headroom to lower limit)
%  A "small" margin (near 0) or negative margin means that bus-phase
%  is congested and the flexible players must act to relieve it.
%% ============================================================

% Confirm the margin interpretation exactly.
margin_upper = vmax - vstar_phys;      % nPhys x 1
margin_lower = vstar_phys - vmin;      % nPhys x 1
assert(norm(b_shF - [margin_upper; margin_lower], inf) < 1e-9, ...
    'b_shF does not match [vmax - vstar_phys; vstar_phys - vmin].');

% Map each physical (non-NaN) non-slack phase slot back to (bus, phase).
physicalFull = idx.nonSlackPhase(valid_ns).';   % full 39-slot physical indices
busOf   = ceil(physicalFull/3);
phaseOf = mod(physicalFull-1,3) + 1;
phChar  = 'abc';
playerBuses = [5 7 12 13];                 % 645, 671, 611, 652

tol_small = 0.02;                          % "near a limit" threshold [pu]

fprintf('\n=== Congestion diagnostics (load_scale = %.3f) ===\n', load_scale);
fprintf('vstar_phys range: [%.4f, %.4f] p.u.  (limits 0.95 / 1.05)\n', ...
    min(vstar_phys), max(vstar_phys));

n_small_up  = nnz(margin_upper < tol_small);
n_small_lo  = nnz(margin_lower < tol_small);
n_neg_up    = nnz(margin_upper < 0);
n_neg_lo    = nnz(margin_lower < 0);
fprintf('Rows with margin < %.3f pu : upper %d, lower %d  (total %d)\n', ...
    tol_small, n_small_up, n_small_lo, n_small_up + n_small_lo);
fprintf('Rows already violated (<0) : upper %d, lower %d\n', ...
    n_neg_up, n_neg_lo);

% List the tight lower-voltage rows (sag toward vmin) — usually the
% binding ones under heavy load — and flag player-adjacent buses.
fprintf('\nTight LOWER-limit rows (vstar_phys near vmin):\n');
[~, ord] = sort(margin_lower, 'ascend');
for r = ord(1:min(15,numel(ord))).'
    tag = '';
    if ismember(busOf(r), playerBuses), tag = '  <-- player bus'; end
    fprintf('  bus %2d phase %c : v = %.4f  margin_lo = %+.4f%s\n', ...
        busOf(r), phChar(phaseOf(r)), vstar_phys(r), margin_lower(r), tag);
end

% Sanity target reminder for the DyNECT / OSQP stage.
fprintf(['\nTarget ~10-30 small-margin rows near player buses [5 7 12 13].\n' ...
         'If too few: raise load_scale or add bus_extra at feeder ends.\n' ...
         'If infeasible (many <0 that players cannot reach): lower load_scale.\n']);

%% Use a path relative to the MATLAB project
export_dir = fullfile(pwd, "data");

if ~exist(export_dir, "dir")
    mkdir(export_dir);
end

export_file = fullfile(export_dir, "network_to_dynect.mat");

% Provenance so the DyNECT / OSQP side knows how this was stressed.
stress_meta = struct( ...
    'load_scale',   load_scale, ...
    'vmin',         vmin(1), ...
    'vmax',         vmax(1), ...
    'vstar_min',    min(vstar_phys), ...
    'vstar_max',    max(vstar_phys), ...
    'm_sh',         m_sh);

save(export_file, ...
    "A_F_game", ...
    "b_shF", ...
    "player_ids", ...
    "n_physical", ...
    "block_start", ...
    "block_end", ...
    "coordinate_convention", ...
    "constraint_convention", ...
    "stress_meta", ...
    "-v7");

fprintf('\nExported network data to:\n%s\n', export_file);

%% comparison
% 
% figure(1)
% 
% phnames = ['a' 'b' 'c'];
% for ph=1:3
% 
%     subplot(3,2,(ph-1)*2+1)
%     plot(1:13, v_testfeeder(rw(1:n,ph)), 'ko', 1:13, v_linearized(rw(1:n,ph)), 'k*');
%     xlim([0 14])
% 	set(gca,'XTick',[1 13])
% 
%     ylabel(sprintf('phase %c',phnames(ph)), 'FontWeight','bold')
% 
% 	if ph==1
% 		title('magnitudes {v_i} [pu]')
% 	end
% 
%     subplot(3,2,(ph-1)*2+2)
%     plot(1:13, t_testfeeder(rw(1:n,ph))/pi*180, 'ko', 1:13, t_linearized(rw(1:n,ph))/pi*180, 'k*')
%     xlim([0 14])
% 	set(gca,'XTick',[1 13])
% 
% 	if ph==1
% 		title('angles \theta_i [deg]')
% 	end
% 
% end

