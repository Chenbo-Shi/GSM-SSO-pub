function res = R_func(p, tau, L, kr)
% R_FUNC Evaluate spherical Bessel/Hankel functions and related Riccati forms.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   res = R_func(p, tau, L, kr)
%
%   This helper function is designed for use in layered-sphere scattering
%   computations. It returns either the spherical Bessel / Hankel function
%   of order L, or the derivative of the corresponding Riccati-Bessel /
%   Riccati-Hankel function divided by kr.
%
%   INPUTS
%     p   = 1 -> spherical Bessel j_L(kr) and its Riccati derivative
%           4 -> spherical Hankel h_L^(2)(kr) and its Riccati derivative
%     tau = 1 -> return spherical function j_L or h_L^(2)
%           2 -> return dR/d(kr), where R = z * f_L(z) and z = kr
%     L   : (integer) order of the spherical function (L = 0,1,2,...)
%     kr  : argument (scalar or array), must be non-zero
%
%   OUTPUT
%     res : same shape as kr, containing the requested quantity

    % Basic input check to avoid division by zero
    if any(kr == 0, 'all')
        error('R_func:ZeroArgument', ...
              'kr must not contain zero to avoid singularities.');
    end

    % Dispatch based on type selector p
    switch p
        case 1  % spherical Bessel j_L
            switch tau
                case 1
                    % Return j_L(kr)
                    res = sph_besselj(L, kr);
                case 2
                    % Return (d/d(kr))[R_j(kr)] / kr, where R_j = kr * j_L(kr)
                    res = dR_j(L, kr) ./ kr;
                otherwise
                    error('R_func:InvalidTau', ...
                          'tau must be 1 (function) or 2 (Riccati derivative).');
            end

        case 4  % spherical Hankel h_L^(2)
            switch tau
                case 1
                    % Return h_L^(2)(kr)
                    res = sph_hankel2(L, kr);
                case 2
                    % Return (d/d(kr))[R_h(kr)] / kr, where R_h = kr * h_L^(2)(kr)
                    res = dR_h(L, kr) ./ kr;
                otherwise
                    error('R_func:InvalidTau', ...
                          'tau must be 1 (function) or 2 (Riccati derivative).');
            end

        otherwise
            error('R_func:InvalidP', ...
                  'p must be 1 (spherical Bessel j_L) or 4 (spherical Hankel h_L^(2)).');
    end

end

%% ===== Local helper functions =====

function j = sph_besselj(v, z)
% SPH_BESSELJ Spherical Bessel function of the first kind j_v(z).
% Uses the relation:
%   j_v(z) = sqrt(pi/(2z)) * J_{v+1/2}(z)
    j = sqrt(pi ./ (2 .* z)) .* besselj(v + 0.5, z);
end

function h2 = sph_hankel2(v, z)
% SPH_HANKEL2 Spherical Hankel function of the second kind h_v^(2)(z).
% Uses the relation:
%   h_v^(2)(z) = sqrt(pi/(2z)) * H_{v+1/2}^{(2)}(z)
    h2 = sqrt(pi ./ (2 .* z)) .* besselh(v + 0.5, 2, z);
end

function dR = dR_j(v, z)
% DR_J Derivative of the Riccati-Bessel function R_j(z) = z * j_v(z).
% Uses the identity:
%   d/dz [z j_v(z)] = (v + 1) * j_v(z) - z * j_{v+1}(z)
    jv  = sph_besselj(v,   z);
    jv1 = sph_besselj(v+1, z);
    dR  = (v + 1) .* jv - z .* jv1;
end

function dR = dR_h(v, z)
% DR_H Derivative of the Riccati-Hankel function R_h(z) = z * h_v^(2)(z).
% Uses the identity:
%   d/dz [z h_v^(2)(z)] = (v + 1) * h_v^(2)(z) - z * h_{v+1}^{(2)}(z)
    hv  = sph_hankel2(v,   z);
    hv1 = sph_hankel2(v+1, z);
    dR  = (v + 1) .* hv - z .* hv1;
end
