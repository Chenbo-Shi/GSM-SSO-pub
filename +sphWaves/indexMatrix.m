function indexMatrix = indexMatrix(maxDegreeL)
% indexMatrix
% Construct the indexing matrix used for ordering spherical-wave modes
% in scattering matrix (S-matrix) formulations.
% (c) 2025, Chenbo Shi, UESTC in China
%
% SYNTAX:
%   indexMatrix = indexMatrix(maxDegreeL)
%
% INPUT:
%   maxDegreeL  – maximum spherical-harmonic degree L used
%
% OUTPUT:
%   indexMatrix – 5 × N matrix encoding mode ordering:
%       row 1 : degree  L
%       row 2 : order   M
%       row 3 : sigma index (1 or 2)
%       row 4 : tau   index (1 or 2)
%       row 5 : alpha ordering index
%
% DESCRIPTION:
%   This function constructs the standard unified indexing for vector
%   spherical wave functions used in electromagnetic scattering.  
%   The construction follows a block pattern:
%       (L, M) pairs → duplicated by sigma and tau (TE/TM patterns)
%   The final alpha-index enforces the canonical ordering defined in
%   analytical spherical-wave expansions.
%
% NOTE:
%   Entries corresponding to (M = 0 AND sigma = 1) are removed because
%   they are overwritten by the alpha-index convention.

%% --- Allocate matrix ---
indexMatrix = nan(5, 2 * maxDegreeL * (maxDegreeL + 3));

%% --- Row 1: degree L for each mode ---
tempL = triu(ones(maxDegreeL + 1, 1) * (0:maxDegreeL));
indexMatrix(1, :) = kron(tempL(tempL ~= 0), ones(4, 1));

%% --- Row 2: order M ---
tempM = (0:maxDegreeL).' * ones(1, maxDegreeL + 1);
tempM = tempM(triu(true(maxDegreeL + 1)));   % keep upper triangular part
indexMatrix(2, :) = kron(tempM(2:end), ones(4, 1));

%% --- Row 3: sigma index (1 1 2 2 repeating) ---
indexMatrix(3, :) = repmat([1 1 2 2], 1, maxDegreeL * (maxDegreeL + 3) / 2);

%% --- Row 4: tau index (1 2 1 2 repeating) ---
indexMatrix(4, :) = repmat([1 2 1 2], 1, maxDegreeL * (maxDegreeL + 3) / 2);

%% --- Row 5: alpha ordering index ---
indexMatrix(5, :) = ...
    2 * ( indexMatrix(1, :).^2 + indexMatrix(1, :) - 1 + ...
          ((-1).^indexMatrix(3, :)) .* indexMatrix(2, :) ) ...
    + indexMatrix(4, :);

%% --- Remove overwritten entries according to alpha-index rules ---
removeIdx = (indexMatrix(2, :) == 0 & indexMatrix(3, :) == 1);
indexMatrix(:, removeIdx) = [];

end
