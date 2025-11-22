function [Ndegree, Lmax] = minLebedevDegree(ka, c1, c2)
% MINLEBEDEVDEGREE
%   Estimate the minimum Lebedev quadrature degree based on electrical size
%   (ka) and user-defined precision constants (c1, c2).
%
% Inputs (optional)
%   c1  ~ tuning constant for estimating the minimum quadrature degree.
%         Default: 2. Higher values increase quadrature precision but also
%         computational cost.
%   c2  ~ similar tuning constant. Default: 3.
%
% Outputs
%   Ndegree ~ recommended Lebedev quadrature degree
%   Lmax    ~ corresponding Lmax order if spherical-wave evaluation and
%             T-matrix decomposition are used
%
% Notes
%   - Numerical formulas follow the referenced paper (see README).
%   - No computational logic has been altered from the original version.
%
% (c) 2025, Chenbo Shi, UESTC in China
% -------------------------------------------------------------------------

nInputs = nargin;

% Substitute default values if c1 or c2 are not provided
if nInputs < 2
    c1 = 2;
end
if nInputs < 3
    c2 = 3;
end

% Use formulas from the paper (see README) to compute Lmax and Ndegree
Lmax    = ceil(ka + c1*ka.^(1/3) + c2);
Ndegree = ceil(4/3 * (Lmax + 1).^2);

end
