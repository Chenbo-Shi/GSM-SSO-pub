function [T, Phi, Psi, Rho] = GenSSO_ana(deg, r, epsi, mu, k0)
% GENSca_Anisophic_optimal Generate scattering matrices for multilayer anisotropic spheres.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [T, PHI, PSI, RHO] = GENSca_Anisophic_optimal(DEG, R, EPSI, MU, K0)
%
%   This function generates the full set of scattering information for
%   a spherically layered anisotropic structure.
%
%   INPUTS:
%     DEG  : maximum degree (angular order) used to build index matrix A
%     R    : 1 x Nlay vector of layer interface radii
%     EPSI : 2 x Nlay array of permittivity parameters
%            EPSI(1,:) -> transverse permittivity (epsilon_t)
%            EPSI(2,:) -> radial     permittivity (epsilon_r)
%     MU   : 2 x Nlay array of permeability parameters
%            MU(1,:)   -> transverse permeability (mu_t)
%            MU(2,:)   -> radial     permeability (mu_r)
%     K0   : wavenumber in vacuum
%
%   OUTPUTS:
%     T    : transition (scattering) matrix (diagonal form)
%     Phi  : inner-to-outer transmission matrix (diagonal form)
%     Psi  : outer-to-inner transmission matrix (diagonal form)
%     Rho  : reflection matrix (diagonal form)
%
%   Internal notes:
%   - Nlay = number of layers.
%   - A    = index matrix mapping (degree/order) to wave index.
%   - Nwav = total number of waves/modes.
%   - For each L, S = 2*(2*L+1) - 1 = 4L + 1 is the total number of
%     degenerate modes (both TE and TM) in the block.

% Number of layers
Nlay = length(r);

% Index matrix (user-supplied function)
A = sphWaves.indexMatrix(deg);
Nwav = size(A, 2);

% Material parameters
et = epsi(1, :);  % transverse permittivity
er = epsi(2, :);  % radial permittivity
mt = mu(1, :);    % transverse permeability
mr = mu(2, :);    % radial permeability

% Precompute layer-dependent quantities (used many times)
Z            = sqrt(mt ./ et);       % impedance-like quantity per layer
ratio_mt_mr  = mt ./ mr;             % mt/mr for effective order of Lg
ratio_et_er  = et ./ er;             % et/er for effective order of Lh
k_layer      = sqrt(et .* mt) * k0;  % effective k in each layer: k = sqrt(et*mt)*k0

%% Forward iteration to compute T matrix
T      = zeros(Nwav, Nlay);
T_last = zeros(Nwav, 1);

alpha = 1;
while alpha <= Nwav
    L = A(1, alpha);
    % Total degeneracy (TE + TM) for this L
    S = 2 * (2 * L + 1) - 1;    % = 4L + 1
    Lfac = L * (L + 1);         % appears repeatedly

    for ilay = 0 : Nlay - 2
        idx   = ilay + 1;       % layer index
        idxp1 = ilay + 2;       % next layer index

        % Radial arguments for this interface at r(idx)
        x_i = k_layer(idx)   * r(idx);
        y_i = k_layer(idxp1) * r(idx);

        % Impedances of current and next layer
        Z_i   = Z(idx);
        Z_ip1 = Z(idxp1);

        % Effective orders for anisotropic medium (L_g and L_h)
        Lg_i   = sqrt(ratio_mt_mr(idx)   * Lfac + 0.25) - 0.5;
        Lh_i   = sqrt(ratio_et_er(idx)   * Lfac + 0.25) - 0.5;
        Lg_ip1 = sqrt(ratio_mt_mr(idxp1) * Lfac + 0.25) - 0.5;
        Lh_ip1 = sqrt(ratio_et_er(idxp1) * Lfac + 0.25) - 0.5;
       
        if ilay == 0
            % First interface: no previous reflection
            F_1i = (1 / Z_i) * utilities.R_func(1, 2, Lg_i, x_i) ./ utilities.R_func(1, 1, Lg_i, x_i);
            F_2i = Z_i       * utilities.R_func(1, 2, Lh_i, x_i) ./ utilities.R_func(1, 1, Lh_i, x_i);
        else
            % Subsequent interfaces: include previous T_last
            num1 = utilities.R_func(1, 2, Lg_i, x_i) + T_last(alpha)   .* utilities.R_func(4, 2, Lg_i, x_i);
            den1 = utilities.R_func(1, 1, Lg_i, x_i) + T_last(alpha)   .* utilities.R_func(4, 1, Lg_i, x_i);
            F_1i = (num1 ./ Z_i) ./ den1;

            num2 = utilities.R_func(1, 2, Lh_i, x_i) + T_last(alpha+1) .* utilities.R_func(4, 2, Lh_i, x_i);
            den2 = utilities.R_func(1, 1, Lh_i, x_i) + T_last(alpha+1) .* utilities.R_func(4, 1, Lh_i, x_i);
            F_2i = (num2 .* Z_i) ./ den2;
        end

        % TE-like block (even indices in this block)
        T(alpha : 2 : alpha + S, ilay+1) = - ...
            ( Z_ip1 * utilities.R_func(1, 1, Lg_ip1, y_i) .* F_1i - utilities.R_func(1, 2, Lg_ip1, y_i) ) ./ ...
            ( Z_ip1 * utilities.R_func(4, 1, Lg_ip1, y_i) .* F_1i - utilities.R_func(4, 2, Lg_ip1, y_i) );

        % TM-like block (odd indices in this block)
        T(alpha+1 : 2 : alpha+1+S, ilay+1) = - ...
            ( (1 / Z_ip1) * utilities.R_func(1, 1, Lh_ip1, y_i) .* F_2i - utilities.R_func(1, 2, Lh_ip1, y_i) ) ./ ...
            ( (1 / Z_ip1) * utilities.R_func(4, 1, Lh_ip1, y_i) .* F_2i - utilities.R_func(4, 2, Lh_ip1, y_i) );

        % Update T_last for the entire degenerate block
        T_last(alpha      : 2 : alpha + S)   = T(alpha,   ilay+1);
        T_last(alpha + 1  : 2 : alpha+1 + S) = T(alpha+1, ilay+1);
    end

    % Move to next L-block
    alpha = alpha + S + 1;
end

%% Backward iteration to compute Psi matrix (outer-to-inner transmission)
Psi      = zeros(Nwav, Nlay);
Psi_last = ones(Nwav, 1);

alpha = 1;
while alpha <= Nwav
    L = A(1, alpha);
    S = 2 * (2 * L + 1) - 1;
    Lfac = L * (L + 1);

    for ilay = Nlay-2 : -1 : 0
        idx   = ilay + 1;
        idxp1 = ilay + 2;

        x_i = k_layer(idx)   * r(idx);
        y_i = k_layer(idxp1) * r(idx);

        Z_i   = Z(idx);
        Z_ip1 = Z(idxp1);

        Lg_i   = sqrt(ratio_mt_mr(idx)   * Lfac + 0.25) - 0.5;
        Lh_i   = sqrt(ratio_et_er(idx)   * Lfac + 0.25) - 0.5;
        Lg_ip1 = sqrt(ratio_mt_mr(idxp1) * Lfac + 0.25) - 0.5;
        Lh_ip1 = sqrt(ratio_et_er(idxp1) * Lfac + 0.25) - 0.5;

        if ilay == 0
            % Innermost interface: denominator does not include T(:, ilay)
            Psi(alpha : 2 : alpha + S, idx) = Psi_last(alpha) .* ...
                ( utilities.R_func(1, 1, Lg_ip1, y_i) + T(alpha,   idx) .* utilities.R_func(4, 1, Lg_ip1, y_i) ) ./ ...
                  utilities.R_func(1, 1, Lg_i, x_i);

            Psi(alpha+1 : 2 : alpha+1+S, idx) = (Z_i / Z_ip1) * Psi_last(alpha+1) .* ...
                ( utilities.R_func(1, 1, Lh_ip1, y_i) + T(alpha+1, idx) .* utilities.R_func(4, 1, Lh_ip1, y_i) ) ./ ...
                  utilities.R_func(1, 1, Lh_i, x_i);
        else
            % Other interfaces: denominator includes T(:, ilay)
            Psi(alpha : 2 : alpha + S, idx) = Psi_last(alpha) .* ...
                ( utilities.R_func(1, 1, Lg_ip1, y_i) + T(alpha,   idx)   .* utilities.R_func(4, 1, Lg_ip1, y_i) ) ./ ...
                ( utilities.R_func(1, 1, Lg_i, x_i)   + T(alpha,   idx)   .* utilities.R_func(4, 1, Lg_i,   x_i) );

            Psi(alpha+1 : 2 : alpha+1+S, idx) = (Z_i / Z_ip1) * Psi_last(alpha+1) .* ...
                ( utilities.R_func(1, 1, Lh_ip1, y_i) + T(alpha+1, idx)   .* utilities.R_func(4, 1, Lh_ip1, y_i) ) ./ ...
                ( utilities.R_func(1, 1, Lh_i, x_i)   + T(alpha+1, ilay) .* utilities.R_func(4, 1, Lh_i,   x_i) );
        end

        Psi_last(alpha      : 2 : alpha + S)   = Psi(alpha,   idx);
        Psi_last(alpha + 1  : 2 : alpha+1 + S) = Psi(alpha+1, idx);
    end

    alpha = alpha + S + 1;
end

%% Backward iteration to compute Rho matrix (reflection)
Rho    = zeros(Nwav, Nlay);
G_last = zeros(Nwav, 1);

alpha = 1;
while alpha <= Nwav
    L = A(1, alpha);
    S = 2 * (2 * L + 1) - 1;
    Lfac = L * (L + 1);

    for ilay = Nlay-2 : -1 : 0
        idx   = ilay + 1;
        idxp1 = ilay + 2;

        x_i = k_layer(idx)   * r(idx);
        y_i = k_layer(idxp1) * r(idx);

        Z_i   = Z(idx);
        Z_ip1 = Z(idxp1);

        Lg_i   = sqrt(ratio_mt_mr(idx)   * Lfac + 0.25) - 0.5;
        Lh_i   = sqrt(ratio_et_er(idx)   * Lfac + 0.25) - 0.5;
        Lg_ip1 = sqrt(ratio_mt_mr(idxp1) * Lfac + 0.25) - 0.5;
        Lh_ip1 = sqrt(ratio_et_er(idxp1) * Lfac + 0.25) - 0.5;

        if ilay == Nlay - 2
            % Outermost interface initialization
            F_1i = (1 / Z_ip1) * utilities.R_func(4, 2, Lg_ip1, y_i) ./ utilities.R_func(4, 1, Lg_ip1, y_i);
            F_2i = Z_ip1       * utilities.R_func(4, 2, Lh_ip1, y_i) ./ utilities.R_func(4, 1, Lh_ip1, y_i);
        else
            % Subsequent interfaces (going inward)
            num1 = G_last(alpha)   .* utilities.R_func(1, 2, Lg_ip1, y_i) + utilities.R_func(4, 2, Lg_ip1, y_i);
            den1 = G_last(alpha)   .* utilities.R_func(1, 1, Lg_ip1, y_i) + utilities.R_func(4, 1, Lg_ip1, y_i);
            F_1i = (num1 ./ Z_ip1) ./ den1;

            num2 = G_last(alpha+1) .* utilities.R_func(1, 2, Lh_ip1, y_i) + utilities.R_func(4, 2, Lh_ip1, y_i);
            den2 = G_last(alpha+1) .* utilities.R_func(1, 1, Lh_ip1, y_i) + utilities.R_func(4, 1, Lh_ip1, y_i);
            F_2i = (num2 .* Z_ip1) ./ den2;
        end

        Rho(alpha : 2 : alpha + S, idx) = - ...
            ( Z_i * utilities.R_func(4, 1, Lg_i, x_i) .* F_1i - utilities.R_func(4, 2, Lg_i, x_i) ) ./ ...
            ( Z_i * utilities.R_func(1, 1, Lg_i, x_i) .* F_1i - utilities.R_func(1, 2, Lg_i, x_i) );

        Rho(alpha+1 : 2 : alpha+1+S, idx) = - ...
            ( (1 / Z_i) * utilities.R_func(4, 1, Lh_i, x_i) .* F_2i - utilities.R_func(4, 2, Lh_i, x_i) ) ./ ...
            ( (1 / Z_i) * utilities.R_func(1, 1, Lh_i, x_i) .* F_2i - utilities.R_func(1, 2, Lh_i, x_i) );

        G_last(alpha      : 2 : alpha + S)   = Rho(alpha,   idx);
        G_last(alpha + 1  : 2 : alpha+1 + S) = Rho(alpha+1, idx);
    end

    alpha = alpha + S + 1;
end

%% Forward iteration to compute Phi matrix (inner-to-outer transmission)
Phi      = zeros(Nwav, Nlay);
Phi_last = ones(Nwav, 1);

alpha = 1;
while alpha <= Nwav
    L = A(1, alpha);
    S = 2 * (2 * L + 1) - 1;
    Lfac = L * (L + 1);

    for ilay = 0 : Nlay - 2
        idx   = ilay + 1;
        idxp1 = ilay + 2;

        x_i = k_layer(idx)   * r(idx);
        y_i = k_layer(idxp1) * r(idx);

        Z_i   = Z(idx);
        Z_ip1 = Z(idxp1);

        Lg_i   = sqrt(ratio_mt_mr(idx)   * Lfac + 0.25) - 0.5;
        Lh_i   = sqrt(ratio_et_er(idx)   * Lfac + 0.25) - 0.5;
        Lg_ip1 = sqrt(ratio_mt_mr(idxp1) * Lfac + 0.25) - 0.5;
        Lh_ip1 = sqrt(ratio_et_er(idxp1) * Lfac + 0.25) - 0.5;

        if ilay == Nlay - 2
            % Outermost interface: denominator uses only outer region
            Phi(alpha : 2 : alpha + S, idx) = Phi_last(alpha) .* ...
                ( Rho(alpha,   idx) .* utilities.R_func(1, 1, Lg_i, x_i) + utilities.R_func(4, 1, Lg_i, x_i) ) ./ ...
                  utilities.R_func(4, 1, Lg_ip1, y_i);

            Phi(alpha+1 : 2 : alpha+1+S, idx) = (Z_ip1 / Z_i) * Phi_last(alpha+1) .* ...
                ( Rho(alpha+1, idx) .* utilities.R_func(1, 1, Lh_i, x_i) + utilities.R_func(4, 1, Lh_i, x_i) ) ./ ...
                  utilities.R_func(4, 1, Lh_ip1, y_i);
        else
            % Interior interfaces: denominator includes Rho(:, ilay+2)
            Phi(alpha : 2 : alpha + S, idx) = Phi_last(alpha) .* ...
                ( Rho(alpha,   idx)   .* utilities.R_func(1, 1, Lg_i,   x_i) + utilities.R_func(4, 1, Lg_i,   x_i) ) ./ ...
                ( Rho(alpha,   idxp1) .* utilities.R_func(1, 1, Lg_ip1, y_i) + utilities.R_func(4, 1, Lg_ip1, y_i) );

            Phi(alpha+1 : 2 : alpha+1+S, idx) = (Z_ip1 / Z_i) * Phi_last(alpha+1) .* ...
                ( Rho(alpha+1, idx)   .* utilities.R_func(1, 1, Lh_i,   x_i) + utilities.R_func(4, 1, Lh_i,   x_i) ) ./ ...
                ( Rho(alpha+1, idxp1) .* utilities.R_func(1, 1, Lh_ip1, y_i) + utilities.R_func(4, 1, Lh_ip1, y_i) );
        end

        Phi_last(alpha      : 2 : alpha + S)   = Phi(alpha,   idx);
        Phi_last(alpha + 1  : 2 : alpha+1 + S) = Phi(alpha+1, idx);
    end

    alpha = alpha + S + 1;
end

%% Extract diagonal (final) matrices from the last/first columns
% These are diagonal matrices in the wave index basis.
T   = diag(T(:,   end-1));  % keep same column as original implementation
Psi = diag(Psi(:, 1));
Phi = diag(Phi(:, end-1));
Rho = diag(Rho(:, 1));

end
