function [Vm, Va, conv] = nrpf(Ybus, Pspec, Qspec, slack, Vslack, Vm0, Va0)
% NRPF  Minimal polar Newton-Raphson AC power flow (all non-slack buses = PQ).
%       Used only to validate the linear input-output map against true physics.
%
%   [Vm, Va, conv] = nrpf(Ybus, Pspec, Qspec, slack, Vslack, Vm0, Va0)
%
%   Pspec, Qspec : n x 1 specified injections (slack entries ignored).
%   slack        : slack bus index, held at Vm=Vslack, Va=0.
%   Vm0, Va0     : optional flat/warm start (default flat 1.0 / 0).
%   conv         : true if converged.

    n = size(Ybus,1);
    if nargin < 6 || isempty(Vm0), Vm0 = ones(n,1); end
    if nargin < 7 || isempty(Va0), Va0 = zeros(n,1); end
    Vm = Vm0; Va = Va0;
    Vm(slack) = Vslack; Va(slack) = 0;

    pq = setdiff(1:n, slack)';           % all non-slack buses are PQ
    G = real(Ybus); B = imag(Ybus);
    tol = 1e-12; conv = false;

    for it = 1:50
        V  = Vm .* exp(1j*Va);
        Sc = V .* conj(Ybus * V);
        dP = Pspec - real(Sc);
        dQ = Qspec - imag(Sc);
        mism = [dP(pq); dQ(pq)];
        if max(abs(mism)) < tol, conv = true; break; end

        % Full polar Jacobian (dense, tiny systems)
        th = Va; Vmag = Vm;
        H = zeros(n); N = zeros(n); M = zeros(n); L = zeros(n);
        for i = 1:n
            for k = 1:n
                a = th(i) - th(k);
                if i ~= k
                    H(i,k) =  Vmag(i)*Vmag(k)*( G(i,k)*sin(a) - B(i,k)*cos(a));
                    N(i,k) =  Vmag(i)*Vmag(k)*( G(i,k)*cos(a) + B(i,k)*sin(a));
                    M(i,k) = -Vmag(i)*Vmag(k)*( G(i,k)*cos(a) + B(i,k)*sin(a));
                    L(i,k) =  Vmag(i)*Vmag(k)*( G(i,k)*sin(a) - B(i,k)*cos(a));
                end
            end
            P_i = real(Sc(i)); Q_i = imag(Sc(i));
            H(i,i) = -Q_i - B(i,i)*Vmag(i)^2;
            N(i,i) =  P_i + G(i,i)*Vmag(i)^2;
            M(i,i) =  P_i - G(i,i)*Vmag(i)^2;
            L(i,i) =  Q_i - B(i,i)*Vmag(i)^2;
        end
        % state = [dTheta(pq); dVmag(pq)/Vmag(pq)] convention -> use standard blocks
        Jpp = H(pq,pq);
        Jpv = N(pq,pq) * diag(1./Vmag(pq));   % dP/dVmag
        Jqp = M(pq,pq);
        Jqv = L(pq,pq) * diag(1./Vmag(pq));   % dQ/dVmag
        J = [Jpp, Jpv; Jqp, Jqv];
        dx = J \ mism;
        npq = numel(pq);
        Va(pq) = Va(pq) + dx(1:npq);
        Vm(pq) = Vm(pq) + dx(npq+1:end);
    end
end
