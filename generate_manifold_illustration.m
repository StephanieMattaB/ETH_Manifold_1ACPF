%% GENERATE_MANIFOLD_ILLUSTRATION
%
% Recreates example_2bus.m's true-manifold-vs-tangent-plane visualization,
% but on a real terminal from our feeder instead of a toy 2-bus system:
% 675.c, the worst (lowest-voltage) terminal at S''. Sweep ITS OWN (p,q)
% around the S'' baseline, holding every other bus-phase's injection fixed,
% and compare:
%   - the TRUE manifold: v(675.c) from an independent nonlinear AC solve
%     at each (p,q) grid point (a genuine re-solve, not an approximation)
%   - the LOCAL TANGENT PLANE: v_pred = v*'' + dv/dp*dp + dv/dq*dq, using
%     S_v's own row/column for 675.c (from the S'' tangent plane already
%     derived in derive_shared_constraint.m -- not re-fit here)
%
% This is the direct visual answer to "why does linear error grow under
% stress": the true manifold curves away from the flat tangent plane, and
% the sweep range here (+-0.06 pu, comparable to the local capability box)
% is exactly the range Sections 6-7 operate in.
%
% Produces: figures/fig9_manifold_675c.png

clear all
close all
clc

run_ieee13_pipeline;

bus_name = '675'; phase_name = 'c';
bus_idx = find(strcmp(bus_labels, bus_name));
phase_idx = 3;   % c
full_idx = rw(bus_idx, phase_idx);
ns_idx = full_idx - 3;   % index into the 36-row non-slack space

idx = sh.idx;
v0 = Spp.ac.v(full_idx);
row = ns_idx;
col_p = ns_idx;
col_q = numel(idx.nonSlackPhase) + ns_idx;
Svp = sh.S_v(row, col_p);
Svq = sh.S_v(row, col_q);

fprintf('Terminal %s.%s: S'''' baseline v = %.6f pu, dv/dp = %.4f, dv/dq = %.4f\n', ...
    bus_name, phase_name, v0, Svp, Svq);

%% Grid sweep: true manifold (independent AC re-solve at every point)

span = 0.06;   % pu, comparable to the local capability box used in Sections 6-7
ngrid = 13;
dp_range = linspace(-span, span, ngrid);
dq_range = linspace(-span, span, ngrid);
[DP, DQ] = meshgrid(dp_range, dq_range);

V_true = zeros(ngrid, ngrid);
V_tangent = v0 + Svp*DP + Svq*DQ;

Vslack = v_slack .* exp(1j*t_slack);
V0_ac = Spp.ac.V;

for i = 1:ngrid
    for j = 1:ngrid
        p_pert = p_S2; q_pert = q_S2;
        p_pert(full_idx) = p_pert(full_idx) + DP(i,j);
        q_pert(full_idx) = q_pert(full_idx) + DQ(i,j);
        Sspec = -(p_pert + 1j*q_pert); Sspec(1:3) = 0;
        ac = solve_acpf_3ph_rect(Y, (1:3)', Vslack, Sspec, V0_ac, ac_opts);
        V_true(i,j) = ac.v(full_idx);
    end
end

err = abs(V_true - V_tangent);
fprintf('Manifold-vs-tangent error over the swept region: max %.4e, RMS %.4e pu\n', ...
    max(err(:)), sqrt(mean(err(:).^2)));

%% Plot: true manifold (mesh) vs tangent plane (transparent surf), anchor marked

figure('visible','off','position',[0,0,1400,600]);

subplot(1,2,1);
mesh(DP, DQ, V_true, 'EdgeColor', 'k'); hold on;
surf(DP, DQ, V_tangent, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'FaceColor', [0.85 0.2 0.2]);
plot3(0, 0, v0, 'k.', 'MarkerSize', 30);
xlabel('\Delta p_{675c} [pu]'); ylabel('\Delta q_{675c} [pu]'); zlabel('|V_{675c}| [pu]');
title('True nonlinear manifold (mesh) vs local tangent plane (red), anchored at S''''');
view(-35, 20); grid on;

subplot(1,2,2);
surf(DP, DQ, err, 'EdgeColor', 'none');
xlabel('\Delta p_{675c} [pu]'); ylabel('\Delta q_{675c} [pu]'); zlabel('|v_{true} - v_{tangent}| [pu]');
title(sprintf('Manifold-vs-tangent error (max %.2e pu over +-%.2f pu sweep)', max(err(:)), span));
colorbar; view(-35, 20); grid on;

print(fullfile('figures','fig9_manifold_675c.png'), '-dpng', '-r150');
fprintf('Saved figures/fig9_manifold_675c.png\n');
