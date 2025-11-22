function [degrees, order, Lmax, Npw, Nsph] = getLebedevDegrees(goal)
% GETLEBEDEVDEGREES
%   Utility for retrieving Lebedev quadrature degrees and the corresponding
%   numbers of plane waves and spherical waves required.
%
% Inputs (optional)
%   goal   ~ target Lebedev quadrature degree. If provided, the closest
%            available degree is selected and only the corresponding
%            single-entry values are returned. If omitted, the full tables
%            (vectors) are returned.
%
% Outputs
%   degrees ~ degree of the Lebedev quadrature
%   order   ~ order of the spherical waves integrated precisely when the
%             corresponding quadrature degree is used
%   Lmax    ~ exact spherical-harmonic order supported
%   Npw     ~ number of plane waves required
%   Nsph    ~ equivalent number of spherical waves (see README for details)
%
% Notes
%   - Numerical values and behavior are kept identical to the original
%     implementation; only documentation and formatting have been updated.
%
% (c) 2025, Chenbo Shi, UESTC in China
%
% -------------------------------------------------------------------------

% (Extracted from bin.getLebedevSphere.m:)
degrees = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, ...
    230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, ...
    1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];

order = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,...
    77,83,89,95,101,107,113,119,125,131];

% Order of spherical waves treated exactly
Lmax = 1:length(order);

% Number of plane waves required
Npw  = 2 * degrees;

% Equivalent number of spherical waves
Nsph = 2 * cumsum(order);

% If "goal" is given, return the closest degree
if nargin > 0
    difference = abs(degrees - goal);
    [~, pos] = min(difference);
    degrees = degrees(pos);
    order   = order(pos);
    Lmax    = Lmax(pos);
    Npw     = Npw(pos);
    Nsph    = Nsph(pos);
end

end
