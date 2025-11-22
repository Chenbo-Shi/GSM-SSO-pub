function [T, Phi, Psi, Rho] = GenSSO_ode_parallel(deg, r, epsi_t, epsi_r, mu_t, mu_r, k0, Reltol)
% GenSca_Anisophic_ode (parallel version)
% Generate scattering matrices for multilayer anisotropic spheres
% using a numerical ODE solver with parallelization over angular blocks.
% The stiffness of the radial ODE system must be assessed by the user.
% For stiff regimes, substitute every ode45 solver call with ode15s.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [T, PHI, PSI, RHO] = GenSca_Anisophic_ode(DEG, R, EPSI_T, EPSI_R, MU_T, MU_R, K0, RELTOL)
%
%   This function generates the full set of scattering information for
%   a spherically layered anisotropic structure by solving radial ODEs.
%   Parallelization is performed over angular-order blocks (L blocks),
%   so each worker independently processes one block and returns its
%   local contribution, which is then assembled on the client.
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

    psiFun  = @(l, x) x .* sbesselj(l, x);
    xiFun   = @(l, x) x .* sbesselh(l, x);
    dpsiFun = @(l, x) (l + 1) .* sbesselj(l, x) - x .* sbesselj(l + 1, x);
    dxiFun  = @(l, x) (l + 1) .* sbesselh(l, x) - x .* sbesselh(l + 1, x);

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

    %% Build angular-order blocks (L-blocks) for parallel processing

    alpha = 1;
    blockStarts = [];
    L_list      = [];
    S_list      = [];

    while alpha <= Nwav
        L = A(1, alpha);
        S = 2 * (2 * L + 1) - 1;  % number of modes in this (TE+TM) block minus 1

        blockStarts(end+1, 1) = alpha; %#ok<AGROW>
        L_list(end+1, 1)      = L;     %#ok<AGROW>
        S_list(end+1, 1)      = S;     %#ok<AGROW>

        alpha = alpha + S + 1;
    end

    Nblocks = numel(L_list);

    % Local block results (each element is a column vector of length S+1)
    T_blocks   = cell(Nblocks, 1);
    Phi_blocks = cell(Nblocks, 1);
    Psi_blocks = cell(Nblocks, 1);
    Rho_blocks = cell(Nblocks, 1);

    %% Parallel loop over blocks

    parfor b = 1:Nblocks
        Lb      = L_list(b);
        Sb      = S_list(b);
        % Each block returns a contiguous vector of length Sb+1
        [T_blocks{b}, Phi_blocks{b}, Psi_blocks{b}, Rho_blocks{b}] = ...
            solveOneDegreeBlock( ...
                Lb, r, et, er, mt, mr, Det, Dmt, ...
                k0, rb, ra, Z0, ...
                psiFun, xiFun, dpsiFun, dxiFun, ...
                Nlay, Sb, options ...
            );
    end

    %% Assemble global vectors from blocks

    T_vec   = zeros(Nwav, 1);
    Phi_vec = zeros(Nwav, 1);
    Psi_vec = zeros(Nwav, 1);
    Rho_vec = zeros(Nwav, 1);

    for b = 1:Nblocks
        alpha0 = blockStarts(b);    % starting index in global vector
        S      = S_list(b);
        idx    = alpha0 : alpha0 + S;

        T_vec(idx)   = T_blocks{b};
        Phi_vec(idx) = Phi_blocks{b};
        Psi_vec(idx) = Psi_blocks{b};
        Rho_vec(idx) = Rho_blocks{b};
    end

    % Return diagonal matrices
    T   = diag(T_vec);
    Psi = diag(Psi_vec);
    Phi = diag(Phi_vec);
    Rho = diag(Rho_vec);
end

% ======================================================================= %
%                           LOCAL HELPER FUNCTIONS                        %
% ======================================================================= %

function [T_block, Phi_block, Psi_block, Rho_block] = solveOneDegreeBlock( ...
    L, r, et, er, mt, mr, Det, Dmt, ...
    k0, rb, ra, Z0, ...
    psiFun, xiFun, dpsiFun, dxiFun, ...
    Nlay, S, options)
% Solve both forward (T, Psi) and backward (Rho, Phi) problems
% for one given angular degree L, and return local block vectors.
%
% The block is of length S+1 = 4*L + 2, containing interleaved TE/TM:
%   indices 1,3,5,...  -> TE-like
%   indices 2,4,6,...  -> TM-like

    blockLen   = S + 1;
    T_block    = zeros(blockLen, 1);
    Phi_block  = zeros(blockLen, 1);
    Psi_block  = zeros(blockLen, 1);
    Rho_block  = zeros(blockLen, 1);

    % Common material data at first and last layers
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

    %% Forward integration: compute T_block (TE/TM) and Psi_block

    % NOTE: The inner loop "for l = 1:L" follows the original implementation.
    % The final values used to build scattering coefficients are based on L.
    for l = 1:L
        % --- TE-like branch (g) forward ---
        for ilay = 1 : Nlay - 2
            mt_i   = mt{ilay};
            et_i   = et{ilay};
            mt_ip  = mt{ilay + 1};
            et_ip  = et{ilay + 1};
            mr_ip  = mr{ilay + 1};
            er_ip  = er{ilay + 1};
            Dmt_ip = Dmt{ilay + 1};
            Det_ip = Det{ilay + 1};

            s1 = @(z) mt_ip(z) ./ mr_ip(z) * (l * (l + 1)) ./ z.^2;
            s2 = @(z) et_ip(z) ./ er_ip(z) * (l * (l + 1)) ./ z.^2;

            p1 = @(z) -Dmt_ip(z) ./ mt_ip(z);
            q1 = @(z) k0^2 * mt_ip(z) .* et_ip(z) - s1(z);

            if ilay == 1
                % Outer boundary condition at rb
                g0        = 1;
                bc_factor = (Z0 / Zb) * dpsiFun(l, x1) / psiFun(l, x1);
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
                g1 = g(1); % store value at inner side of first interval
            end

            % --- TM-like branch (h) forward ---
            p2 = @(z) -Det_ip(z) ./ et_ip(z);
            q2 = @(z) k0^2 * mt_ip(z) .* et_ip(z) - s2(z);

            if ilay == 1
                h0        = 1;
                bc_factor = (Zb / Z0) * dpsiFun(l, x1) / psiFun(l, x1);
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
                dh1 = dh(1) / et_ip(rb);
            end
        end

        % Use g at ra to build outer scattering (TE-like)
        factor_te = mt_f(0) / kf / mt_ip(ra) * dg(end) / g(end);
        % TE indices inside this block: 1,3,5,...,blockLen-1
        T_block(1:2:blockLen) = ...
            -(factor_te * psiFun(L, y2) - dpsiFun(L, y2)) ./ ...
              (factor_te * xiFun(L,  y2) - dxiFun(L,  y2));

        Psi_block(1:2:blockLen) = g1 / g(end) * ra / rb .* ...
            (utilities.R_func(1, 1, L, y2) + T_block(1) * utilities.R_func(4, 1, L, y2)) ./ ...
             utilities.R_func(1, 1, L, x1);

        % Use h at ra to build outer scattering (TM-like)
        factor_tm = et_f(0) / kf / et_ip(ra) * dh(end) / h(end);
        % TM indices inside this block: 2,4,6,...,blockLen
        T_block(2:2:blockLen) = ...
            -(factor_tm * psiFun(L, y2) - dpsiFun(L, y2)) ./ ...
              (factor_tm * xiFun(L,  y2) - dxiFun(L,  y2));

        Psi_block(2:2:blockLen) = dh1 / dh(end) * ra / rb .* et_ip(ra) .* ...
            (utilities.R_func(1, 2, L, y2) + T_block(2) * utilities.R_func(4, 2, L, y2)) ./ ...
             utilities.R_func(1, 2, L, x1);
    end

    %% Backward integration: compute Rho_block and Phi_block

    for l = 1:L
        % --- TE-like branch (g) backward ---
        for ilay = Nlay - 1 : -1 : 2
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
                % Outer boundary at ra (backward)
                g0        = 1;
                bc_factor = (Z0 / Zf) * dxiFun(l, y2) / xiFun(l, y2);
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

            % --- TM-like branch (h) backward ---
            p2 = @(z) -Det_i(z) ./ et_i(z);
            q2 = @(z) k0^2 * mt_i(z) .* et_i(z) - s2(z);

            if ilay == Nlay - 1
                h0        = 1;
                bc_factor = (Zf / Z0) * dxiFun(l, y2) / xiFun(l, y2);
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

        % Use g at rb to build inner scattering (TE-like)
        factor_te = mt_b(0) / kb / mt_i(rb) * dg(end) / g(end);
        Rho_block(1:2:blockLen) = ...
            -(factor_te * xiFun(L, x1) - dxiFun(L, x1)) ./ ...
              (factor_te * psiFun(L, x1) - dpsiFun(L, x1));

        Phi_block(1:2:blockLen) = g1 / g(end) * rb / ra .* ...
            (Rho_block(1) * utilities.R_func(1, 1, L, x1) + utilities.R_func(4, 1, L, x1)) ./ ...
             utilities.R_func(4, 1, L, y2);

        % Use h at rb to build inner scattering (TM-like)
        factor_tm = et_b(0) / kb / et_i(rb) * dh(end) / h(end);
        Rho_block(2:2:blockLen) = ...
            -(factor_tm * xiFun(L, x1) - dxiFun(L, x1)) ./ ...
              (factor_tm * psiFun(L, x1) - dpsiFun(L, x1));

        Phi_block(2:2:blockLen) = dh1 / dh(end) * rb / ra .* et_i(rb) .* ...
            (Rho_block(2) * utilities.R_func(1, 2, L, x1) + utilities.R_func(4, 2, L, x1)) ./ ...
             utilities.R_func(4, 2, L, y2);
    end
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
