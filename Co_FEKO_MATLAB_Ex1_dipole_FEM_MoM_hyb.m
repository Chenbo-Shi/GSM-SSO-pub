%% Frequency sweep: dipole FEM–MoM hybrid example and field post-processing
% (c) 2025, Chenbo Shi, UESTC in China
%
%   This script computes the generalized scattering matrix of a dipole
%   antenna enclosed inside a cylindrical dielectric cavity, using a
%   hybrid FEM–MoM formulation.
%
%   Note:
%   This script does *not* evaluate the behavior of the same structure
%   when embedded in a dielectric sphere. That analysis is performed
%   separately in script Ex5.

modelname = 'dipole_FEM_MoM_hyb';

% Parse basic cards from existing FEKO .pre model
[SF, IN, DP, EG, PS, CG, AW, DI, FP] = ...
    feko.parse_pre_fields(['model/' modelname '.pre']);

%% Build media and port descriptions
MediumTable = feko.buildDItable(DI);  % Material definitions table

% DP is an N×1 string array of port definitions
% buildPortTables also snaps very small numbers to zero
portTables = feko.buildPortTables(DP);
keysList   = keys(portTables);

for i = 1:numel(keysList)
    portName = keysList{i};
    fprintf('--- Port: %s ---\n', portName);
    disp(portTables(portName));
end

% Annotate port information using AW card(s)
portInfo = feko.annotatePortsFromAW(portTables, AW);

% Get model scale from SF card
[scale, unit] = feko.getScaleFromSF(SF); %#ok<ASGLU>

%% Frequency definition
freqs = 1.5e9 : 0.1e9 : 6e9;
a     = 36 * scale;            % Radius of minimum enclosing sphere
k0    = 2 * pi * freqs / 3e8;
ka    = k0 * a;
Z0    = 376.73031;

nFreqs = length(freqs);

%% Build header cards for new .pre files
IN = ['IN   8 1055  "' ['model/' modelname '.cfm"']];
% Example of using direct sparse solver (commented out by default):
% CG = 'CG: 21 :  : -1 :  :  :  :  :  :  :  :  : 0';

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

%% Frequency sweep loop
for nf = 1:nFreqs
    f0 = freqs(nf);

    % Estimate Lebedev degree and Lmax from electrical size
    [nPW_est, Lmax] = lebedev.minLebedevDegree(ka(nf));
    % nPW_est = 120;  % user-defined override
    % Lmax    = 17;   % user-defined override

    nDegree_est = lebedev.getLebedevDegrees(nPW_est);
    nDegree     = nDegree_est;

    %% Upload Lebedev quadrature
    quadrature = lebedev.getLebedevSphere(nDegree);

    % Theta, phi for Lebedev quadrature points
    [~, th, ph] = lebedev.cart2sph(quadrature.x, quadrature.y, quadrature.z);

    r2d = @(x) 180 * x / pi;  % Radians to degrees

    % Plane-wave excitation directions
    % Rotate angles to match the definition used for scattering dyadics
    thPW = r2d(pi - [th; th]);
    phPW = r2d(pi + [ph; ph]);
    pol  = r2d([pi*ones(nDegree,1); -pi/2*ones(nDegree,1)]);  % [theta phi]
    PW   = [thPW, phPW, pol];

    % Far-field observation points
    FF      = r2d([th ph]);
    FF_xyz  = [quadrature.x, quadrature.y, quadrature.z]; %#ok<NASGU>
    pwInfo  = struct('PW', PW, 'FF', FF);

    %% Count total number of propagating port modes at this frequency
    % (used to size matrices and S-parameter file names)
    nFarfields = nDegree;        %#ok<NASGU>
    nWaves     = 2 * nDegree;

    portInfo = feko.countPropagatingModesAtFmax( ...
        portInfo, MediumTable, scale, f0);

    nPortmodes = 0;
    for i = 1:length(portInfo)
        nPortmodes = nPortmodes + size(portInfo(i).modesList, 1);
    end

    % Write a new FEKO .pre file for this frequency and solution setup
    feko.write_prefeko([preFEKOfileName '.pre'], cards, f0, portInfo, pwInfo);

    % Solve in FEKO (may be memory intensive with multiple RHS)
    dos(['runfeko ' preFEKOfileName '.pre --execute-prefeko -np all']);
    clc;

    %% Load S-parameters from .snp file
    Gam = sparameters( ...
        [preFEKOfileName '_Spara.s' num2str(nPortmodes) 'p']).Parameters;

    % Input power scaling, for lossless systems:
    % Pin_scale = 1 - |Gamma|^2 per port.
    % (For lossy systems, the absorbed power must be subtracted separately.)
    Pin_scale = diag(eye(size(Gam)) - Gam' * Gam);

    % Allocate scattering dyadic matrices
    thisRadi = nan(nWaves, nPortmodes, 1);
    thisSdya = nan(nWaves, nWaves,     1);

    %% Load radiated fields from FEKO (.ffe files)
    thisRadi(:, 1, :) = feko.readFEKORadiatingFields( ...
        preFEKOfileName, nDegree);

    parfor n = 2:nPortmodes
        thismode = ['(' num2str(n-1) ')'];
        thisRadi(:, n, :) = feko.readFEKORadiatingFields( ...
            preFEKOfileName, nDegree, thismode);
    end

    Radi = thisRadi / sqrt(Z0);

    %% Load scattered fields from FEKO (.ffe files)
    thisSdya(:, 1, :) = feko.readFEKOScattertingFields( ...
        preFEKOfileName, nDegree);

    parfor n = 2:nWaves
        thisPW = ['(' num2str(n-1) ')'];
        thisSdya(:, n, :) = feko.readFEKOScattertingFields( ...
            preFEKOfileName, nDegree, thisPW);
    end

    % Scale factor -j k0 / (4π) for scattering dyadic
    Sdyadic = thisSdya .* reshape(-1j * k0(nf) / (4*pi), 1, 1, []);

    %% Vector spherical harmonics on Lebedev points (frequency independent)
    indexMatrix      = sphWaves.indexMatrix(Lmax);
    [P, Nsphw]       = feko.getcomplexSphHarmonic(indexMatrix, th, ph);
    w                = diag(quadrature.w);  % Lebedev weights (diagonal matrix)

    % Build T-matrix from radiated fields
    P_th       = squeeze(P(:,:,2));
    P_ph       = squeeze(P(:,:,3));
    P_lebedev  = [P_th P_ph];  % Complex vector spherical harmonics

    % Orthogonality check (not used later, but kept for completeness)
    D = P_th * w * P_th' + P_ph * w * P_ph'; %#ok<NASGU>

    T = conj(P_lebedev) * blkdiag(w, w) * Radi * diag(sqrt(Pin_scale));
    U = conj(P_lebedev) * blkdiag(w, w) * Sdyadic * blkdiag(w, w) * P_lebedev.';

    % Scattering matrix in spherical wave basis
    S  = eye(Nsphw) + 2 * U;
    Sg = [Gam T.'; T S];  % Generalized scattering matrix

    %% Save data for this frequency
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
