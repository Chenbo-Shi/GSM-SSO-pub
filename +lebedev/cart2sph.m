function [r, theta, phi] = cart2sph(x, y, z)
% (c) 2025, Chenbo Shi, UESTC in China
% CART2SPH  Transform Cartesian coordinates (x,y,z) to spherical (r,theta,phi).
%
%   [r, theta, phi] = CART2SPH(x, y, z)
%
% Inputs
%   x : x-coordinate, numeric array (same size as y,z)
%   y : y-coordinate, numeric array
%   z : z-coordinate, numeric array
%
% Outputs
%   r     : radial coordinate, same size as input
%   theta : polar angle (0 <= theta <= pi)
%   phi   : azimuthal angle (-pi < phi <= pi)
%
% Notes
%   - Vectorized implementation; assumes x, y, z have identical size.
%   - Numerical formulas match the original behavior exactly.

% Input consistency check
if ~isequal(size(x), size(y), size(z))
    error('cart2sph:SizeMismatch', ...
          'Inputs x, y, z must have identical size.');
end

% Cartesian → spherical
r     = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./ r);
phi   = atan2(y, x);
end
