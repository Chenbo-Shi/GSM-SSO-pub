%% Homogeneous, anisotropic, lossy medium example (Sec. IV-B)
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script computes the generalized scattering matrix of a 5-mode horn
% antenna embedded in a concentric, homogeneous, anisotropic, lossy
% spherical structure. The horn is modeled by its generalized scattering
% matrix in free space, which must be precomputed and stored in
% preserved_data/CoSim_Ex0_Horn_Liang_Script.
%
% For each frequency, the script:
%   1) Loads the precomputed generalized scattering matrix Sg of the horn.
%   2) Partitions Sg into port and spherical-wave sub-blocks.
%   3) Uses analytical spherical shell operators (GenSSO_ana) for an
%      anisotropic medium to obtain the modified reflection matrix
%      \tilde{\Gamma}.
%   4) Stores selected diagonal and coupling entries of the modified
%      reflection matrix versus frequency.
%
% The resulting curves correspond to Sec. IV-B of the paper.

%% Problem definition
freqs = 3.2 : 0.05 : 3.8;     % Frequency range [GHz]
nfreq = numel(freqs);
c     = 299792458;            % Speed of light [m/s]

% Radii of concentric spherical layers [m] (from inner to outer).
% The last "inf" represents free space.
r = [150 180 inf] * 1e-3;

% Relative permittivity and permeability:
%   epsi_t, mu_t: transverse components
%   epsi_r, mu_r: radial components
epsi_t = [1 5 1];
epsi_r = [1 2 1];
mu_t   = [1 3 1];
mu_r   = [1 1 1];

% Output array for selected elements of the modified reflection matrix
% Columns: [Γ11 Γ22 Γ33 Γ44 Γ55 Γ45]
Spara_new = zeros(nfreq, 6);

dataDir  = fullfile('preserved_data', 'CoSim_Ex0_Horn_Liang_Script');
baseName = 'Horn_Liang_Script_data';

%% Frequency sweep
for i = 1:nfreq
    f = freqs(i)  % Frequency in GHz (for reference)

    % Load precomputed generalized scattering matrix Sg
    dataFile = fullfile(dataDir, sprintf('%s(%d).mat', baseName, i));
    load(dataFile, 'Sg', 'nPortmodes', 'Lmax'); %#ok<LOAD>

    % Size of spherical-wave sub-block
    Nswf = 2 * Lmax * (Lmax + 2);

    % Partition Sg into sub-blocks:
    %   Sg = [ Gam  R
    %          T    S ]
    Gam = Sg(1:nPortmodes,               1:nPortmodes);              % Port-to-port
    R   = Sg(1:nPortmodes,               nPortmodes+1:end);          % Port-to-SWF
    T   = Sg(nPortmodes+1:end,           1:nPortmodes);              % SWF-to-port
    U   = (Sg(nPortmodes+1:end, nPortmodes+1:end) - eye(Nswf)) / 2;  % Internal SWF block

    % Analytical spherical shell operator for anisotropic, lossy medium
    omega = 2 * pi * f * 1e9;  % Angular frequency [rad/s]
    [t, Phi, Psi, Rho] = utilities.GenSSO_ana( ... %#ok<ASGLU>
        Lmax, r, [epsi_t; epsi_r], [mu_t; mu_r], omega / c);

    % Modified reflection matrix (ports only)
    Gam_new = Gam + R/2 * Rho / (eye(Nswf) - U * Rho) * T;

    % Store selected entries as a row vector
    Spara_new(i, :) = [ ...
        Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
        Gam_new(4,4), Gam_new(5,5), Gam_new(4,5)];
end

%% Plot results: 20*log10|Gamma_tilde|
fig = figure; %#ok<NASGU>
set(gcf, 'Color', 'w');

hold on;
h1 = plot(freqs, 20*log10(abs(Spara_new)), 'LineWidth', 1);

ax = gca;
set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8, ...
    'Box', 'on', ...
    'LineWidth', 0.5, ...
    'GridLineStyle', ':', ...
    'YMinorTick', 'off', ...
    'MinorGridLineStyle', 'none', ...
    'GridColor', 'k');

xlim([freqs(1) freqs(end)]);
xticks(freqs);

ylabel('$20\log_{10}|\tilde{\mathbf{\Gamma}}|$', ...
    'Interpreter', 'latex');
xlabel('Frequency (GHz)');

grid on;

legend(ax, h1, ...
    {'$\tilde{\mathbf{\Gamma}}_{11}$', '$\tilde{\mathbf{\Gamma}}_{22}$', ...
     '$\tilde{\mathbf{\Gamma}}_{33}$', '$\tilde{\mathbf{\Gamma}}_{44}$', ...
     '$\tilde{\mathbf{\Gamma}}_{55}$', '$\tilde{\mathbf{\Gamma}}_{45}$'}, ...
    'NumColumns', 3, ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8, ...
    'Box', 'on', ...
    'Location', 'southwest');

set(gca, 'GridLineStyle', ':', 'GridColor', 'k');
