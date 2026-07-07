function A_vd = build_Avd(Ybus, V0, A0)
% BUILD_AVD  Bolognani-Dorfler tangent matrix A_vd of the power-flow manifold.
%
%   A_vd = build_Avd(Ybus, V0, A0)
%
%   Returns the 2n-by-2n matrix A_vd such that, in deviations from the
%   linearization point (V0, A0),
%
%         A_vd * [dVmag; dTheta] = [dP; dQ].
%
%   This is exactly the left block of Agrid = [A_vd  -I] that example_case14.m
%   assembles inline (S. Bolognani, F. Dorfler, Allerton 2015). It is factored
%   out here so the GNE input-output map can reuse it.
%
%   DIFFERENCE FROM example_case14.m: that script hardcodes
%   PP = Rmatrix(ones,zeros) because it linearizes at a FLAT voltage profile.
%   Here PP = Rmatrix(V0, A0), which is the correct polar->rectangular Jacobian
%   at a GENERAL operating point. Pass the true AC solution as (V0, A0) to get
%   a map that is first-order exact at the actual operating point.
%
%   Ordering convention (matches example_case14.m):
%     rows  1..n   -> dP ,   rows  n+1..2n -> dQ
%     cols  1..n   -> dVmag, cols  n+1..2n -> dTheta

    n = numel(V0);

    U    = V0 .* exp(1j*A0);        % complex voltage phasor at lin. point
    J0   = Ybus * U;               % complex current injection
    UU   = bracket(diag(U));       % 2n x 2n
    JJ   = bracket(diag(conj(J0)));% 2n x 2n
    NN   = Nmatrix(2*n);           % [I 0; 0 -I]
    YY   = bracket(Ybus);          % 2n x 2n
    PP   = Rmatrix(V0, A0);        % polar->rect Jacobian at (V0,A0)

    A_vd = (JJ + UU*NN*YY) * PP;   % 2n x 2n
end
