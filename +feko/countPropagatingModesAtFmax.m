function portInfo = countPropagatingModesAtFmax(portInfo, MediumTable, scale, freq_max, varargin)
% countPropagatingModesAtFmax
%   Determine propagating modes at the maximum frequency for each port.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   portInfo = countPropagatingModesAtFmax(portInfo, MediumTable, scale, freq_max, ...)
%
%   For each portInfo(i), this function fills:
%       portInfo(i).modesList : each row = [typeCode, m, n, polDeg]
%           typeCode : 1 = TE, 2 = TM, 3 = TEM
%           m, n     : mode indices
%                      (for coax TEM, [3, 0, 0, 0])
%           polDeg   : polarization flag (0 or 90 degrees)
%
%   Only truly propagating modes are included (i.e., kz not purely imaginary).
%
%   INPUTS
%     portInfo   : struct array, each element at least has fields
%                  .typeCode  (1=rectangular, 2=circular, 3=coaxial)
%                  .mediumLabel
%                  .XYZ       (3×N point coordinates per port)
%     MediumTable: table with columns Label, EpsRel, Sigma, TanD
%     scale      : geometric scale factor for XYZ
%     freq_max   : maximum frequency of interest (Hz)
%
%   OPTIONAL NAME–VALUE PAIRS
%     'mMax'     : maximum m index (default 10, rectangular only)
%     'nMax'     : maximum n index (default 10, rectangular only)
%
%   NOTES
%     - Rectangular waveguide (typeCode = 1):
%         Rectangular TE/TM modes up to (mMax, nMax) are scanned and
%         propagating ones are retained.
%     - Circular waveguide (typeCode = 2):
%         Not supported in this simplified version; no modes are added.
%     - Coaxial line (typeCode = 3):
%         Only fundamental TEM mode is considered; if a valid inner and
%         outer conductor exist (r_out > r_in > 0), one TEM mode [3,0,0,0]
%         is added.
%


% ---------------- options ----------------
p = inputParser;
p.addParameter('mMax', 10, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('nMax', 10, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.parse(varargin{:});
mMax = p.Results.mMax;
nMax = p.Results.nMax;

% ---------------- constants ----------------
epsilon0  = 8.85e-12;
omega_max = 2 * pi * freq_max;
k0_max    = omega_max / 3e8;

absTol = 1e-9;
relTol = 1e-9;

for i = 1:numel(portInfo)
    type = portInfo(i).typeCode;

    % --- material parameters from MediumTable (default: free space) ---
    idx = strcmp(MediumTable.Label, portInfo(i).mediumLabel);
    if any(idx)
        kRow = find(idx, 1, 'first');
        eps_r = MediumTable.EpsRel(kRow);
        sigma = MediumTable.Sigma(kRow);
        tand  = MediumTable.TanD(kRow);
    else
        eps_r = 1;
        sigma = 0;
        tand  = 0;
    end

    % complex permittivity and wavenumber in the filling material
    eps_c = eps_r * (1 - 1i * tand) - 1i * sigma / (epsilon0 * omega_max);
    kg    = k0_max * sqrt(eps_c);

    modes = zeros(0, 4);   % [typeCode, m, n, polDeg]

    switch type
        case 1  % ===== rectangular waveguide =====
            XYZ = portInfo(i).XYZ * scale;
            S1  = XYZ(:, 1);
            S2  = XYZ(:, 2);
            S3  = XYZ(:, 3);
            a   = norm(S2 - S1);  % x dimension
            b   = norm(S3 - S1);  % y dimension

            for n = 0:nMax
                for m = 0:mMax
                    if m + n == 0
                        continue;  % skip (0,0)
                    end
                    kc2 = (m * pi / a)^2 + (n * pi / b)^2;
                    kz  = sqrt(kg^2 - kc2);
                    if isNonPureImag(kz, absTol, relTol)
                        % When both m and n are nonzero, TE/TM are distinct.
                        % When m==0 or n==0, only one polarization exists.
                        if (m == 0) || (n == 0)
                            % Only TE or TM is physically present; here we
                            % keep TE (type=1) as representative.
                            modes(end+1, :) = [1, m, n, 0]; %#ok<AGROW>
                        else
                            modes(end+1, :) = [1, m, n, 0]; %#ok<AGROW> % TE
                            modes(end+1, :) = [2, m, n, 0]; %#ok<AGROW> % TM
                        end
                    end
                end
            end

        case 2  % ===== circular waveguide =====
            % Not supported in this simplified version: leave modes empty
            % (no modes added).

        case 3  % ===== coaxial line (TEM only) =====
            XYZ = portInfo(i).XYZ * scale;

            % Convention: S1=center, S2=inner-conductor point, S4=outer-conductor point
            % (the actual ordering of points must match the upstream definition)
            S1 = XYZ(:, 1);
            S2 = XYZ(:, 2);
            S4 = XYZ(:, 4);

            a = norm(S2 - S1);
            b = norm(S4 - S1);

            r_out = max(a, b);
            r_in  = min(a, b);

            % TEM exists if there are two conductors, r_out > r_in > 0
            if (r_out > r_in) && (r_in > 0)
                modes(end+1, :) = [3, 0, 0, 0];  %#ok<AGROW>  % TEM
            end

        otherwise
            % Unknown type: no modes added
    end

    portInfo(i).modesList = modes;
end
end

% ============================================================
% Criterion: kz is considered propagating if it is not purely
% imaginary within given tolerances.
% ============================================================
function tf = isNonPureImag(kz, absTol, relTol)
tf = abs(real(kz)) > max(absTol, relTol * abs(kz));
end
