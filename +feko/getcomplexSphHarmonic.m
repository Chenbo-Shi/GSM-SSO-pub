function [P, Nsphw] = getcomplexSphHarmonic(indexMatrix, th, ph)
% getcomplexSphHarmonic Compute complex vector spherical harmonics.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [P, Nsphw] = getcomplexSphHarmonic(indexMatrix, th, ph)
%
%   This function constructs complex vector spherical harmonics based on a
%   given index matrix and angular coordinates (theta, phi), using the
%   real-valued vector spherical harmonics returned by sphWaves.functionA.
%
%   INPUTS:
%     indexMatrix : 5 x Nsphw index matrix that encodes the ordering and
%                   type of spherical wave functions:
%                     row 1 -> degree l
%                     row 2 -> order  m
%                     row 3 -> sigma index (even/odd)
%                     row 4 -> tau   index (1/2 type)
%                     row 5 -> alpha ordering index
%     th          : column vector of polar angles theta  [M x 1], in radians
%     ph          : column vector of azimuthal angles phi [M x 1], in radians
%
%   OUTPUTS:
%     P     : complex vector spherical harmonics, size [Nsphw x M x 3]
%     Nsphw : number of spherical wave functions (columns of indexMatrix)

    % Number of spherical wave functions
    Nsphw = size(indexMatrix, 2);

    % Real-valued vector spherical harmonics A1, A2 (A3 is not needed here)
    [A1, A2, ~] = sphWaves.functionA( ...
        indexMatrix(1, :).', indexMatrix(2, :).', th, ph);

    % Preallocate real-valued vector spherical harmonics container
    A = nan(size(A1), 'like', A1);

    % Logical masks for tau / sigma combinations
    ind_1e = (indexMatrix(4, :) == 1) & (indexMatrix(3, :) == 2); % tau = 1, sigma = e
    ind_1o = (indexMatrix(4, :) == 1) & (indexMatrix(3, :) == 1); % tau = 1, sigma = o
    ind_2e = (indexMatrix(4, :) == 2) & (indexMatrix(3, :) == 2); % tau = 2, sigma = e
    ind_2o = (indexMatrix(4, :) == 2) & (indexMatrix(3, :) == 1); % tau = 2, sigma = o

    % Build real-valued vector spherical harmonics from A1, A2
    A(ind_1e, :, :) = real(A1(ind_1e, :, :));
    A(ind_1o, :, :) = imag(A1(ind_1o, :, :));
    A(ind_2e, :, :) = real(A2(ind_2e, :, :));
    A(ind_2o, :, :) = imag(A2(ind_2o, :, :));

    % Complex phase factor: -(-j)^(tau - l)
    pre_script = -(-1j) .^ (indexMatrix(4, :).' - indexMatrix(1, :).');
    pre_script = reshape(pre_script, Nsphw, 1, 1);  % for implicit expansion

    % Sort in alpha sequence (row 5 of indexMatrix)
    [~, alpha_seq] = sort(indexMatrix(5, :).');
    A = A(alpha_seq, :, :);

    % Apply complex phase factor to obtain complex vector spherical harmonics
    P = A .* pre_script;
end
