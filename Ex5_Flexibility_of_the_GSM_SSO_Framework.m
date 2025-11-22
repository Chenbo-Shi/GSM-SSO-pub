%% Embedded dipole in spherical dielectric shell: generalized reflection coefficient
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script evaluates the performance of a dipole antenna, already
% encapsulated in a cylindrical dielectric cavity (FEM–MoM hybrid model),
% when it is further embedded in a concentric spherical dielectric shell.
%
% The generalized scattering matrix Sg of the FEM–MoM dipole is obtained
% by the FEKO–MATLAB co-simulation script:
%   Co_FEKO_MATLAB_Ex1_dipole_FEM_MoM_hyb.m
%
% Precomputed data are stored under:
%   preserved_data/CoSim_Ex1_dipole_FEM_MoM_hyb
%
% For each frequency this script:
%   1) Loads Sg, nPortmodes and Lmax.
%   2) Builds the analytical spherical-shell operator Rho.
%   3) Updates the port reflection matrix Γ (including the spherical shell).
%   4) Stores Γ11 and plots its magnitude in dB versus frequency.

%% Frequency sweep [GHz]
fbeg  = 1.5;
fend  = 6.0;
fstep = 0.1;
freqs = fbeg:fstep:fend;
nFreq = numel(freqs);

%% Spherical dielectric shell definition
% Radii of concentric spherical layers [m] (inner radius, outer radius, free space)
r  = [45 60 inf] * 1e-3;
% Relative permittivity and permeability of each region
er = [1 5*(1-0.1j) 1];
mr = [1 1          1];

%% Data location (precomputed generalized scattering matrices)
dataDir  = fullfile('preserved_data', 'CoSim_Ex1_dipole_FEM_MoM_hyb');
baseName = 'dipole_FEM_MoM_hyb_data';

% Storage for updated S-parameter (reflection coefficient of the dipole port)
Spara_new = complex(zeros(nFreq, 1));

%% Frequency-by-frequency update of the reflection coefficient
for i = 1:nFreq
    f = freqs(i)  % Frequency in GHz, kept for clarity

    % Load precomputed Sg, nPortmodes and Lmax
    dataFile = fullfile(dataDir, sprintf('%s(%d).mat', baseName, i));
    load(dataFile, 'Sg', 'nPortmodes', 'Lmax'); %#ok<LOAD>

    % Size of spherical-wave block
    Nfwav = 2 * Lmax * (Lmax + 2);

    % Truncate Sg to the port + spherical-wave subspace (safety)
    Sg = Sg(1:nPortmodes+Nfwav, 1:nPortmodes+Nfwav);

    % Partition Sg:
    %   Sg = [ Gam  R
    %          T    S ]
    Gam = Sg(1:nPortmodes,               1:nPortmodes);              % Port-to-port
    R   = Sg(1:nPortmodes,               nPortmodes+1:end);          % Port-to-SWF
    T   = Sg(nPortmodes+1:end,           1:nPortmodes);              % SWF-to-port
    U   = (Sg(nPortmodes+1:end, nPortmodes+1:end) - eye(Nfwav)) / 2; % Internal SWF block

    % Analytical spherical-shell operator for homogeneous, isotropic, lossy shell
    [~, ~, ~, Rho] = utilities.GenSSO_ana(Lmax, r, [er; er], [mr; mr], 2*pi*f*1e9/3e8);

    % Updated port reflection matrix including the spherical shell
    Gam_new = Gam + R/2 * Rho * inv(eye(Nfwav) - U * Rho) * T;

    % Store reflection coefficient of the (single) dipole port
    Spara_new(i) = Gam_new(1,1);
end

%% Plot |Γ| (dB) versus frequency
figure;
hold on;

h1 = plot(freqs, 20*log10(abs(Spara_new)), 'LineWidth', 1);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'LineWidth', 0.5, ...
    'GridLineStyle', ':', 'yMinorTick', 'off', ...
    'MinorGridLineStyle', 'none', 'GridColor', 'k');

xlim([fbeg fend]);
xticks(1.5:0.5:6);

ylabel('$\Gamma$ (dB)', 'Interpreter', 'latex');
xlabel('Frequency (GHz)');
grid on;

legend(h1, {'FEM + MoM'}, 'NumColumns', 1, ...
   'FontName', 'Times New Roman', 'FontSize', 8, ...
   'Box', 'on', 'Location', 'southeast');

set(gcf, 'Color', 'w');
set(gca, 'GridLineStyle', ':', 'GridColor', 'k');
