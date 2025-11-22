function [A1, A2, A3] = functionA(degreeL, orderM, theta, phi)
% functionA  Vector spherical harmonics A (three orthogonal vector fields)
% (c) 2025, Chenbo Shi, UESTC in China
%
%  INPUTS
%   degreeL : vector of degrees l,           double [N x 1]
%   orderM  : vector of orders  m,           double [N x 1]
%   theta   : polar angle coordinates (rad), double [M x 1]
%   phi     : azimuthal angle coordinates (rad), double [M x 1]
%
%  OUTPUTS
%   A1 : vector spherical harmonic A1, complex double [N x M x 3]
%   A2 : vector spherical harmonic A2, complex double [N x M x 3]
%   A3 : vector spherical harmonic A3, complex double [N x M x 3]
%
%  SYNTAX
%   [A1, A2, A3] = functionA(degreeL, orderM, theta, phi)
%
%  NOTES
%   - This implementation assumes that degreeL and orderM are column
%     vectors of the same length.
%   - Implicit expansion (R2016b or later) is used for broadcasting.

    % --- Basic input checks (lightweight, but helps catch shape issues) ---
    if numel(degreeL) ~= numel(orderM)
        error('functionA:SizeMismatch', ...
              'degreeL and orderM must have the same number of elements.');
    end

    % Ensure column vectors
    degreeL = degreeL(:);
    orderM  = orderM(:);
    theta   = theta(:);
    phi     = phi(:);

    % --- Dimensions ---
    nRows = numel(degreeL);   % number of (l, m) pairs
    nCols = numel(theta);     % number of angular sample points

    % --- Preallocate outputs (N x M x 3) ---
    A1 = complex(zeros(nRows, nCols, 3));
    A2 = complex(zeros(nRows, nCols, 3));
    A3 = complex(zeros(nRows, nCols, 3));

    % --- Precompute constants not depending on theta/phi ---
    % Remove (-1)^m factor (here effectively set to (+1)^m = 1^m)
    constant1 = (1).^orderM ./ sqrt(2 * pi * degreeL .* (degreeL + 1));
    nonZeroM  = (orderM ~= 0);
    constant1(nonZeroM) = sqrt(2) * constant1(nonZeroM);

    constant2 = (1).^orderM ./ sqrt(2 * pi);
    constant2(nonZeroM) = sqrt(2) * constant2(nonZeroM);

    % exp(i m phi), using implicit expansion: [N x 1] * [1 x M] -> [N x M]
    expPart = exp(1i * orderM * phi.');

    % --- Legendre-related quantities ---
    cosTheta = cos(theta);
    [normLegendreSing, normLegendreSingDer, normLegendre] = ...
        sphWaves.legendreComponents(degreeL, orderM, cosTheta);

    % --- Assemble vector spherical harmonics ---
    % A1: purely angular (no radial component)
    %   A1_theta component  -> index 2
    %   A1_phi   component  -> index 3
    A1(:, :, 2) = constant1 .* 1i .* normLegendreSing    .* expPart;
    A1(:, :, 3) = -constant1       .* normLegendreSingDer .* expPart;

    % A2: orthogonal angular combination
    A2(:, :, 2) = -A1(:, :, 3);
    A2(:, :, 3) =  A1(:, :, 2);

    % A3: purely radial part
    A3(:, :, 1) = constant2 .* normLegendre .* expPart;
end
