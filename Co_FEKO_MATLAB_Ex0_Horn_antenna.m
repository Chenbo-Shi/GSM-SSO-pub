%% Frequency-by-frequency simulation using FEKO (loop sweep)
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script:
%   1) Reads a reference FEKO .pre file and parses its cards.
%   2) Builds material and port tables.
%   3) For each frequency:
%       - Selects an appropriate Lebedev quadrature.
%       - Sets up plane-wave excitations and far-field observation points.
%       - Counts propagating port modes and writes a new .pre file.
%       - Runs FEKO, reads S-parameters, radiated fields and scattering fields.
%       - Constructs the dyadic scattering matrix and generalized scattering matrix.
%       - Stores results in a MAT file.
%   4) Optionally deletes temporary FEKO files at the end.
%
%   This script is used to compute the free-space generalized scattering
%   matrix of a horn antenna with five operating modes. This model is used
%   in the "Numerical Validation" section.


%% Load base model and parse .pre cards
modelname = 'Horn_Liang_Script';
[SF, IN, DP, EG, PS, CG, AW, DI, FP] = ...
    feko.parse_pre_fields(['model/' modelname '.pre']);

%% Build material and port tables
MediumTable = feko.buildDItable(DI);  % Material definitions from DI

% DP is an N×1 string array
portTables = feko.buildPortTables(DP);   % Tiny numerical values are thresholded to zero by default
keysList   = keys(portTables);
for i = 1:numel(keysList)
    portName = keysList{i};
    fprintf('--- Port: %s ---\n', portName);
    disp(portTables(portName));
end

portInfo        = feko.annotatePortsFromAW(portTables, AW);
[scale, unit]   = feko.getScaleFromSF(SF); %#ok<ASGLU>

%% Define frequency sweep
freqs = 3.2e9 : 0.05e9 : 3.8e9;
a     = 146.7 * scale;             % Minimum enclosing-sphere radius
k0    = 2 * pi * freqs / 3e8;
ka    = k0 * a;
Z0    = 376.73031;

nFreqs = numel(freqs);

%% Prepare cards for new .pre files
IN = ['IN   8 1055  "' ['model/' modelname '.cfm"']];  % keep original content

cards = struct( ...
    'SF', SF, ...
    'IN', IN, ...
    'DP', DP, ...
    'EG', EG, ...
    'PS', PS, ...
    'CG', CG, ...
    'DI', DI, ...
    'FP', FP);

preFEKOfileName = 'GSM_Sdyadic';

%% Frequency loop: run FEKO and post-process results
for nf = 1:nFreqs
    f0 = freqs(nf);

    % Estimate Lebedev degree and Lmax from electrical size
    [nPW_est, Lmax] = lebedev.minLebedevDegree(ka(nf));
    % nPW_est = 120; % user-defined override (if desired)
    % Lmax    = 17;  % user-defined override (if desired)

    % Map required number of plane waves to Lebedev degree
    nDegree_est = lebedev.getLebedevDegrees(nPW_est);
    nDegree     = nDegree_est;

    %% Lebedev quadrature on the unit sphere
    quadrature = lebedev.getLebedevSphere(nDegree);

    % Spherical coordinates (in radians) of the quadrature points
    [~, th, ph] = lebedev.cart2sph(quadrature.x, quadrature.y, quadrature.z);

    % Radians to degrees conversion
    r2d = @(x) 180 * x / pi;

    %% Plane-wave excitations and far-field observation directions
    % Excitation via plane waves (two polarizations per direction)
    thPW = r2d(pi - [th; th]);                    % Rotate to our definition of scattering dyadics
    phPW = r2d(pi + [ph; ph]);
    pol  = r2d([pi*ones(nDegree,1); -pi/2*ones(nDegree,1)]); % [theta, phi] polarization
    PW   = [thPW, phPW, pol];

    % Far-field observation directions
    FF    = r2d([th ph]);
    FF_xyz = [quadrature.x, quadrature.y, quadrature.z]; %#ok<NASGU> % (not used downstream, kept for potential diagnostics)
    pwInfo = struct( ...
        'PW', PW, ...
        'FF', FF);

    %% Count propagating port modes at this frequency
    nFarfields = nDegree;          %#ok<NASGU> % kept for completeness
    nWaves     = 2 * nDegree;      % Total number of plane waves

    portInfo   = feko.countPropagatingModesAtFmax(portInfo, MediumTable, scale, f0);
    nPortmodes = 0;
    for i = 1:numel(portInfo)
        nPortmodes = nPortmodes + size(portInfo(i).modesList, 1);
    end

    %% Write FEKO .pre file for this frequency and run FEKO
    feko.write_prefeko([preFEKOfileName '.pre'], cards, f0, portInfo, pwInfo);

    % Multiple RHS may lead to memory issues in some cases
    dos(['runfeko ' preFEKOfileName '.pre --execute-prefeko -np all']);
    clc;

    %% Load S-parameters from .snp file
    Gam = sparameters([preFEKOfileName '_Spara.s' num2str(nPortmodes) 'p']).Parameters;
    % For a lossless system: Pin_scale = diag(I - Gamma^H Gamma)
    Pin_scale = diag(eye(size(Gam)) - Gam' * Gam);

    % Allocation for dyadic scattering quantities
    thisRadi = nan(nWaves, nPortmodes, 1);
    thisSdya = nan(nWaves, nWaves, 1);

    %% Load radiating far fields (port excitations)
    thisRadi(:, 1, :) = feko.readFEKORadiatingFields(preFEKOfileName, nDegree);
    parfor n = 2:nPortmodes
        thismode = ['(' num2str(n-1) ')'];
        thisRadi(:, n, :) = feko.readFEKORadiatingFields( ...
            preFEKOfileName, nDegree, thismode);
    end
    Radi = thisRadi / sqrt(Z0);

    %% Load scattering far fields (plane-wave excitations)
    thisSdya(:, 1, :) = feko.readFEKOScattertingFields(preFEKOfileName, nDegree);
    parfor n = 2:nWaves
        thisPW = ['(' num2str(n-1) ')'];
        thisSdya(:, n, :) = feko.readFEKOScattertingFields( ...
            preFEKOfileName, nDegree, thisPW);
    end

    % Apply scattering dyadic scaling: -j k / (4*pi)
    Sdyadic = thisSdya .* reshape(-1j * k0(nf) / (4*pi), 1, 1, []);

    %% Complex vector spherical harmonics at Lebedev points (frequency-independent)
    indexMatrix     = sphWaves.indexMatrix(Lmax);
    [P, Nsphw]      = feko.getcomplexSphHarmonic(indexMatrix, th, ph);
    w               = diag(quadrature.w);          % Diagonal matrix of Lebedev weights

    % Decompose P into theta and phi components
    P_th        = squeeze(P(:,:,2));
    P_ph        = squeeze(P(:,:,3));
    P_lebedev   = [P_th P_ph];                     % Complex vector spherical harmonics

    % Check orthogonality (optional diagnostic)
    D = P_th * w * P_th' + P_ph * w * P_ph';       %#ok<NASGU>

    %% Build T-matrix and generalized scattering matrix
    % Emission matrix (port → spherical waves)
    T = conj(P_lebedev) * blkdiag(w, w) * Radi * diag(sqrt(Pin_scale));

    % Scattering part (waves → waves) for a single frequency
    U = conj(P_lebedev) * blkdiag(w, w) * Sdyadic * blkdiag(w, w) * P_lebedev.';

    % S (spherical-wave scattering matrix)
    S = eye(Nsphw) + 2 * U;

    % Generalized scattering matrix:
    %   Sg = [Gamma  T^T;
    %         T      S  ]
    Sg = [Gam T.'; T S];

    %% Save results for this frequency
    tempFileName = sprintf('%s_data(%d).mat', modelname, nf);
    save(tempFileName, ...
        'f0', 'Sdyadic', 'Radi', 'Sg', 'nPortmodes', 'Lmax');
end

%% Delete auxiliary FEKO files (optional)
deleteAuxFiles = true;
if deleteAuxFiles
    feko.cleanupFEKOfarfieldFiles(preFEKOfileName);
    feko.cleanupFEKOsnpFiles(preFEKOfileName);
    feko.cleanupFEKOsolverFiles(preFEKOfileName, {'.pre'});
    feko.cleanupFEKOconfigfile();
end
