%% Radially Piecewise Homogeneous Medium (Sec. IV-C)
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script:
%   1) Loads precomputed generalized scattering matrices Sg for a
%      5-mode horn antenna in free space (stored in preserved_data).
%   2) Embeds the antenna into radially piecewise homogeneous spherical
%      shells (isotropic or anisotropic).
%   3) Computes the modified port reflection matrix Γ̃ for each case using
%      an analytical spherical-shell operator.
%   4) Plots 20·log10|Γ̃_ij| for selected ports versus frequency.
%
% The reference free-space generalized scattering matrix Sg is generated
% by the FEKO–MATLAB co-simulation script:
%   Co_FEKO_MATLAB_Ex0_Horn_antenna.m

%% Common settings
freqs = 3.2 : 0.05 : 3.8;     % Frequency range [GHz]
nfreq = numel(freqs);
c     = 299792458;            % Speed of light [m/s]

% Directory and file name pattern for precomputed Sg data
dataDir  = fullfile('preserved_data', 'CoSim_Ex0_Horn_Liang_Script');
baseName = 'Horn_Liang_Script_data';

%% Case 1: Two-layer isotropic shell
% Radii of concentric spherical layers [m] (from inner to outer).
% The last "inf" represents free space.
r_case1      = [150 165 180 inf] * 1e-3;
epsi_t_case1 = [1 4.4 - 0.396j 10 1];
epsi_r_case1 = [1 4.4 - 0.396j 10 1];
mu_t_case1   = [1 1          1  1];
mu_r_case1   = [1 1          1  1];

% Output array for selected elements of the modified reflection matrix
% Columns: [Γ11 Γ22 Γ33 Γ44 Γ55 Γ45]
Spara_case1 = [];

% Frequency sweep
for i = 1:nfreq
    f = freqs(i) % Frequency in GHz (for reference)

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

    % Analytical spherical-shell operator for isotropic, lossy medium
    [~, Phi, Psi, Rho] = utilities.GenSSO_ana( ... %#ok<ASGLU>
        Lmax, r_case1, [epsi_t_case1; epsi_r_case1], ...
        [mu_t_case1; mu_r_case1], 2 * pi * freqs(i) * 1e9 / c);

    % Modified port reflection matrix
    Gam_new = Gam + R / 2 * Rho / (eye(Nswf) - U * Rho) * T;

    % Store selected elements
    Spara_case1 = [Spara_case1; ...
        [Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
         Gam_new(4,4), Gam_new(5,5), Gam_new(4,5)]]; %#ok<AGROW>
end

% Plot results for case 1
fig1 = figure; %#ok<NASGU>
h1   = plot(freqs, 20 * log10(abs(Spara_case1)), 'LineWidth', 1);
hold on;

set(gca, ...
    'FontName',        'Times New Roman', ...
    'FontSize',        8, ...
    'Box',             'on', ...
    'LineWidth',       0.5, ...
    'GridLineStyle',   ':', ...
    'YMinorTick',      'off', ...
    'MinorGridLineStyle', 'none', ...
    'GridColor',       'k');

xlim([freqs(1), freqs(end)]);
xticks(freqs);

ylabel('$20\log_{10}|\tilde{\mathbf{\Gamma}}|$', ...
    'Interpreter', 'latex');
xlabel('Frequency (GHz)');

grid on;

legend(h1, ...
    {'$\tilde{\mathbf{\Gamma}}_{11}$', ...
     '$\tilde{\mathbf{\Gamma}}_{22}$', ...
     '$\tilde{\mathbf{\Gamma}}_{33}$', ...
     '$\tilde{\mathbf{\Gamma}}_{44}$', ...
     '$\tilde{\mathbf{\Gamma}}_{55}$', ...
     '$\tilde{\mathbf{\Gamma}}_{45}$'}, ...
    'NumColumns', 3, ...
    'Interpreter', 'latex', ...
    'FontName',   'Times New Roman', ...
    'FontSize',   8, ...
    'Box',        'on', ...
    'Location',   'southwest');

set(gcf, 'Color', 'w');
set(gca, 'GridLineStyle', ':', 'GridColor', 'k');

%% Case 2: Two-layer anisotropic shell
r_case2      = [150 165 180 inf] * 1e-3;
epsi_t_case2 = [1 4.4 8 1];
epsi_r_case2 = [1 2   1 1];
mu_t_case2   = [1 2.2 5 1];
mu_r_case2   = [1 2.2 2 1];

Spara_case2 = [];

for i = 1:nfreq
    f = freqs(i)

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

    % Analytical spherical-shell operator for anisotropic shell
    [~, Phi, Psi, Rho] = utilities.GenSSO_ana( ... %#ok<ASGLU>
        Lmax, r_case2, [epsi_t_case2; epsi_r_case2], ...
        [mu_t_case2; mu_r_case2], 2 * pi * freqs(i) * 1e9 / 3e8);

    % Modified port reflection matrix
    Gam_new = Gam + R / 2 * Rho / (eye(Nswf) - U * Rho) * T;

    % Store selected elements
    Spara_case2 = [Spara_case2; ...
        [Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
         Gam_new(4,4), Gam_new(5,5), Gam_new(4,5)]]; %#ok<AGROW>
end

% Plot results for case 2
fig2 = figure; %#ok<NASGU>
h2   = plot(freqs, 20 * log10(abs(Spara_case2)), 'LineWidth', 1);
hold on;

set(gca, ...
    'FontName',        'Times New Roman', ...
    'FontSize',        8, ...
    'Box',             'on', ...
    'LineWidth',       0.5, ...
    'GridLineStyle',   ':', ...
    'YMinorTick',      'off', ...
    'MinorGridLineStyle', 'none', ...
    'GridColor',       'k');

xlim([freqs(1), freqs(end)]);
xticks(freqs);

ylabel('$20\log_{10}|\tilde{\mathbf{\Gamma}}|$', ...
    'Interpreter', 'latex');
xlabel('Frequency (GHz)');

grid on;

legend(h2, ...
    {'$\tilde{\mathbf{\Gamma}}_{11}$', ...
     '$\tilde{\mathbf{\Gamma}}_{22}$', ...
     '$\tilde{\mathbf{\Gamma}}_{33}$', ...
     '$\tilde{\mathbf{\Gamma}}_{44}$', ...
     '$\tilde{\mathbf{\Gamma}}_{55}$', ...
     '$\tilde{\mathbf{\Gamma}}_{45}$'}, ...
    'NumColumns', 3, ...
    'Interpreter', 'latex', ...
    'FontName',   'Times New Roman', ...
    'FontSize',   8, ...
    'Box',        'on', ...
    'Location',   'southwest');

set(gcf, 'Color', 'w');
set(gca, 'GridLineStyle', ':', 'GridColor', 'k');
