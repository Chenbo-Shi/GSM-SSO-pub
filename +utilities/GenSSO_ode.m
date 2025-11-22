function [T, Phi, Psi, Rho] = GenSSO_ode(deg, r, epsi_t, epsi_r, mu_t, mu_r, k0, Reltol)
% GenSca_Anisophic_ode
% Generate scattering matrices for multilayer anisotropic spheres
% using a numerical ODE solver.
% The stiffness of the radial ODE system must be assessed by the user.
% For stiff regimes, substitute every ode45 solver call with ode15s.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [T, PHI, PSI, RHO] = GenSca_Anisophic_ode(DEG, R, EPSI_T, EPSI_R, MU_T, MU_R, K0, RELTOL)
%
%   This function generates the full set of scattering information for
%   a spherically layered anisotropic structure by solving radial ODEs.
%
%   INPUTS:
%     DEG     : maximum degree (angular order) used to build index matrix A
%     R       : 1 x Nlay vector of layer interface radii
%     EPSI_T  : function handle or array describing transverse permittivity (epsilon_t)
%     EPSI_R  : function handle or array describing radial     permittivity (epsilon_r)
%     MU_T    : function handle or array describing transverse permeability (mu_t)
%     MU_R    : function handle or array describing radial     permeability (mu_r)
%     K0      : wavenumber in vacuum
%     RELTOL  : relative tolerance of the ODE solver (optional, default 1e-6)
%
%   OUTPUTS:
%     T       : transition (scattering) matrix (diagonal)
%     Phi     : inner-to-outer transmission matrix (diagonal)
%     Psi     : outer-to-inner transmission matrix (diagonal)
%     Rho     : reflection matrix (diagonal)

    if nargin < 8
        Reltol = 1e-6;
    end

    % Spherical Bessel/Hankel and Riccati functions
    sbesselj = @(l, x) sqrt(pi ./ (2 .* x)) .* besselj(l + 0.5, x);
    sbesselh = @(l, x) sqrt(pi ./ (2 .* x)) .* besselh(l + 0.5, 2, x);

    psi  = @(l, x) x .* sbesselj(l, x);
    xi   = @(l, x) x .* sbesselh(l, x);
    dpsi = @(l, x) (l + 1) .* sbesselj(l, x) - x .* sbesselj(l + 1, x);
    dxi  = @(l, x) (l + 1) .* sbesselh(l, x) - x .* sbesselh(l + 1, x);

    % Number of layers and modes
    Nlay = numel(r);
    A    = sphWaves.indexMatrix(deg);
    Nwav = size(A, 2);

    % Derivatives of transverse material parameters (with respect to radius)
    Dmu_t   = derivative(mu_t);
    Depsi_t = derivative(epsi_t);

    Z0 = 1;
    rb = r(1);
    ra = r(Nlay - 1);

    % Convert material descriptors into layer-resolved function handles
    [et, er, mt, mr, Det, Dmt] = deal(cell(1, Nlay));
    for i = 1:Nlay
        et{i}  = get(epsi_t, i);
        er{i}  = get(epsi_r, i);
        mt{i}  = get(mu_t,   i);
        mr{i}  = get(mu_r,   i);
        Dmt{i} = get(Dmu_t,  i);
        Det{i} = get(Depsi_t,i);
    end

    options = odeset('RelTol', Reltol);

    %% Forward integration to obtain T and Psi
    [T_vec, Psi_vec, Rho_vec, Phi_vec] = deal(zeros(Nwav, 1));

    alpha = 1;
    while alpha <= Nwav
        L = A(1, alpha);
        S = 2 * (2 * L + 1) - 1;

        mt_b = mt{1};
        et_b = et{1};
        mt_f = mt{Nlay};
        et_f = et{Nlay};

        Zb = sqrt(mt_b(0) / et_b(0));
        kb = k0 * sqrt(mt_b(0) * et_b(0));
        kf = k0 * sqrt(mt_f(0) * et_f(0));

        x1 = kb * rb;
        y2 = kf * ra;

        for l = 1:L
            % ----- TE-like branch (g) -----
            for ilay = 1 : Nlay - 2
                mt_i    = mt{ilay};
                et_i    = et{ilay};
                mt_ip   = mt{ilay + 1};
                et_ip   = et{ilay + 1};
                mr_ip   = mr{ilay + 1};
                er_ip   = er{ilay + 1};
                Dmt_ip  = Dmt{ilay + 1};
                Det_ip  = Det{ilay + 1};

                s1 = @(z) mt_ip(z) ./ mr_ip(z) * (l * (l + 1)) ./ z.^2;
                s2 = @(z) et_ip(z) ./ er_ip(z) * (l * (l + 1)) ./ z.^2;

                p1 = @(z) -Dmt_ip(z) ./ mt_ip(z);
                q1 = @(z) k0^2 * mt_ip(z) .* et_ip(z) - s1(z);

                % Solve for g in [r(ilay), r(ilay+1)]
                if ilay == 1
                    % Outer boundary condition at rb
                    g0        = 1;
                    bc_factor = (Z0 / Zb) * dpsi(l, x1) / psi(l, x1);
                    dg0       = g0 * k0 * mt_ip(rb) * bc_factor;
                else
                    % Continuity at internal interface
                    g0  = g(end);
                    dg0 = dg(end) * mt_ip(r(ilay)) / mt_i(r(ilay));
                end

                odefun1 = @(rr, y) [y(2);
                                    -p1(rr) .* y(2) - q1(rr) .* y(1)];
                r_span1 = [r(ilay), r(ilay + 1)];
                y01     = [g0; dg0];

                [~, y_sol] = ode45(odefun1, r_span1, y01, options);
                g  = y_sol(:, 1);
                dg = y_sol(:, 2);

                if ilay == 1
                    g1 = g(1);  % store value at inner boundary of first interval
                end

                % ----- TM-like branch (h) -----
                p2 = @(z) -Det_ip(z) ./ et_ip(z);
                q2 = @(z) k0^2 * mt_ip(z) .* et_ip(z) - s2(z);

                if ilay == 1
                    h0        = 1;
                    bc_factor = (Zb / Z0) * dpsi(l, x1) / psi(l, x1);
                    dh0       = h0 * k0 * et_ip(rb) * bc_factor;
                else
                    h0  = h(end);
                    dh0 = dh(end) * et_ip(r(ilay)) / et_i(r(ilay));
                end

                odefun2 = @(rr, y) [y(2);
                                    -p2(rr) .* y(2) - q2(rr) .* y(1)];
                y02     = [h0; dh0];
                [~, y_sol] = ode45(odefun2, r_span1, y02, options);
                h  = y_sol(:, 1);
                dh = y_sol(:, 2);

                if ilay == 1
                    dh1 = dh(1) / et_ip(rb);  % store derivative-related value at rb
                end
            end

            % Use g at ra to build outer scattering (TE)
            factor_te = mt_f(0) / kf / mt_ip(ra) * dg(end) / g(end);
            T_vec(alpha:2:alpha+S) = ...
                -(factor_te * psi(L, y2) - dpsi(L, y2)) ./ ...
                  (factor_te * xi(L,  y2) - dxi(L,  y2));

            Psi_vec(alpha:2:alpha+S) = g1 / g(end) * ra / rb .* ...
                (utilities.R_func(1,1,L,y2) + T_vec(alpha) * utilities.R_func(4,1,L,y2)) ./ ...
                 utilities.R_func(1,1,L,x1);

            % Use h at ra to build outer scattering (TM)
            factor_tm = et_f(0) / kf / et_ip(ra) * dh(end) / h(end);
            T_vec(alpha+1:2:alpha+1+S) = ...
                -(factor_tm * psi(L, y2) - dpsi(L, y2)) ./ ...
                  (factor_tm * xi(L,  y2) - dxi(L,  y2));

            Psi_vec(alpha+1:2:alpha+1+S) = dh1 / dh(end) * ra / rb .* et_ip(ra) .* ...
                (utilities.R_func(1,2,L,y2) + T_vec(alpha+1) * utilities.R_func(4,2,L,y2)) ./ ...
                 utilities.R_func(1,2,L,x1);
        end

        alpha = alpha + S + 1;
    end

    %% Backward integration to obtain Rho and Phi

    alpha = 1;
    while alpha <= Nwav
        L = A(1, alpha);
        S = 2 * (2 * L + 1) - 1;

        mt_b = mt{1};
        et_b = et{1};
        mt_f = mt{Nlay};
        et_f = et{Nlay};

        Zb = sqrt(mt_b(0) / et_b(0));
        Zf = sqrt(mt_f(0) / et_f(0));
        kb = k0 * sqrt(mt_b(0) * et_b(0));
        kf = k0 * sqrt(mt_f(0) * et_f(0));

        x1 = kb * rb;
        y2 = kf * ra;

        for l = 1:L
            % ----- TE-like branch (g, backward) -----
            for ilay = Nlay-1 : -1 : 2
                mt_i   = mt{ilay};
                et_i   = et{ilay};
                mr_i   = mr{ilay};
                er_i   = er{ilay};
                Dmt_i  = Dmt{ilay};
                Det_i  = Det{ilay};
                mt_ip  = mt{ilay + 1};
                et_ip  = et{ilay + 1};

                s1 = @(z) mt_i(z) ./ mr_i(z) * (l * (l + 1)) ./ z.^2;
                s2 = @(z) et_i(z) ./ er_i(z) * (l * (l + 1)) ./ z.^2;

                p1 = @(z) -Dmt_i(z) ./ mt_i(z);
                q1 = @(z) k0^2 * mt_i(z) .* et_i(z) - s1(z);

                if ilay == Nlay - 1
                    % Outer boundary at ra
                    g0        = 1;
                    bc_factor = (Z0 / Zf) * dxi(l, y2) / xi(l, y2);
                    dg0       = g0 * k0 * mt_i(ra) * bc_factor;
                else
                    g0  = g(end);
                    dg0 = dg(end) * mt_i(r(ilay)) / mt_ip(r(ilay));
                end

                odefun1 = @(rr, y) [y(2);
                                    -p1(rr) .* y(2) - q1(rr) .* y(1)];
                r_span1 = [r(ilay), r(ilay - 1)];
                y01     = [g0; dg0];
                [~, y_sol] = ode45(odefun1, r_span1, y01, options);
                g  = y_sol(:, 1);
                dg = y_sol(:, 2);

                if ilay == Nlay - 1
                    g1 = g(1);
                end

                % ----- TM-like branch (h, backward) -----
                p2 = @(z) -Det_i(z) ./ et_i(z);
                q2 = @(z) k0^2 * mt_i(z) .* et_i(z) - s2(z);

                if ilay == Nlay - 1
                    h0        = 1;
                    bc_factor = (Zf / Z0) * dxi(l, y2) / xi(l, y2);
                    dh0       = h0 * k0 * et_i(ra) * bc_factor;
                else
                    h0  = h(end);
                    dh0 = dh(end) * et_i(r(ilay)) / et_ip(r(ilay));
                end

                odefun2 = @(rr, y) [y(2);
                                    -p2(rr) .* y(2) - q2(rr) .* y(1)];
                y02     = [h0; dh0];
                [~, y_sol] = ode45(odefun2, r_span1, y02, options);
                h  = y_sol(:, 1);
                dh = y_sol(:, 2);

                if ilay == Nlay - 1
                    dh1 = dh(1) / et_i(ra);
                end
            end

            % Use g at rb to build inner scattering (TE)
            factor_te = mt_b(0) / kb / mt_i(rb) * dg(end) / g(end);
            Rho_vec(alpha:2:alpha+S) = ...
                -(factor_te * xi(L, x1) - dxi(L, x1)) ./ ...
                  (factor_te * psi(L, x1) - dpsi(L, x1));

            Phi_vec(alpha:2:alpha+S) = g1 / g(end) * rb / ra .* ...
                (Rho_vec(alpha) * utilities.R_func(1,1,L,x1) + utilities.R_func(4,1,L,x1)) ./ ...
                 utilities.R_func(4,1,L,y2);

            % Use h at rb to build inner scattering (TM)
            factor_tm = et_b(0) / kb / et_i(rb) * dh(end) / h(end);
            Rho_vec(alpha+1:2:alpha+1+S) = ...
                -(factor_tm * xi(L, x1) - dxi(L, x1)) ./ ...
                  (factor_tm * psi(L, x1) - dpsi(L, x1));

            Phi_vec(alpha+1:2:alpha+1+S) = dh1 / dh(end) * rb / ra .* et_i(rb) .* ...
                (Rho_vec(alpha+1) * utilities.R_func(1,2,L,x1) + utilities.R_func(4,2,L,x1)) ./ ...
                 utilities.R_func(4,2,L,y2);
        end

        alpha = alpha + S + 1;
    end

    % Return diagonal matrices
    T   = diag(T_vec);
    Psi = diag(Psi_vec);
    Phi = diag(Phi_vec);
    Rho = diag(Rho_vec);
end

function f = get(fhandle, i)
% Wrap a function handle returning vectors into a scalar-per-layer function.
%   INPUT:  fhandle(z) -> vector of length Nlay
%   OUTPUT: f(z)      -> scalar, i-th component
    f = @(z) subsref(fhandle(z), substruct('()', {i}));
end

function df = derivative(fhandle)
%DERIVATIVE Compute symbolic derivative of a single-variable anonymous function.
%
%   df = DERIVATIVE(fhandle)
%   fhandle : anonymous function of one variable, e.g., @(z) z.^2 + sin(z)
%   df      : anonymous function representing the derivative of fhandle.

    syms z
    try
        f_sym = fhandle(z);
    catch
        error(['Input anonymous function must accept a single symbolic ', ...
               'variable and be evaluable at that variable.']);
    end

    df_sym = diff(f_sym, z);
    df     = matlabFunction(df_sym, 'Vars', z);
end
