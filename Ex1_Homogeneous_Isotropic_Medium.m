%% Homogeneous, isotropic, lossy medium example (Sec. IV-A)
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script:
%   - Loads pre-computed generalized scattering matrices Sg for a 5-mode
%     horn antenna in free space (stored in the "preserved_data" folder).
%   - Embeds the antenna into a concentric, layered, homogeneous, isotropic,
%     lossy spherical medium.
%   - Computes the modified port scattering matrix \tilde{Γ} for each
%     frequency.
%   - Produces a single, publication-ready figure of 20*log10(|\tilde{Γ}|).
%
% NOTE:
%   - If the generalized scattering matrices Sg have not been computed yet,
%     please run the script "Co_FEKO_MATLAB_Ex0_Horn_antenna.m" first.

%% Problem definition
freqs = 3.2 : 0.05 : 3.8;     % Frequency range [GHz]
nfreq = numel(freqs);
c     = 299792458;            % Speed of light [m/s]

% Radii of concentric spherical layers [m] (from inner to outer).
% The last layer is free space (represented by "inf").
r    = [150 180 inf] * 1e-3;                 % [m]
epsi = [1 5 * (1 - 0.1j) 1];                 % Relative permittivity
mu   = [1 1 1];                              % Relative permeability

%% Allocate storage for port scattering results
% Spara_new(i,:) = [Γ_11 Γ_22 Γ_33 Γ_44 Γ_55 Γ_45] at frequency freqs(i)
Spara_new = nan(nfreq, 6);

%% Frequency loop: load Sg, embed into layered medium, compute Γ~
for i = 1:nfreq
    f = freqs(i) % Frequency in GHz (kept for clarity if needed later)

    % Load precomputed generalized scattering matrix Sg
    dataFile = fullfile('preserved_data', ...
        'CoSim_Ex0_Horn_Liang_Script', ...
        sprintf('Horn_Liang_Script_data(%d).mat', i));
    load(dataFile, 'Sg', 'nPortmodes', 'Lmax'); %#ok<LOAD>

    % Size of spherical-wave sub-block
    Nswf = 2 * Lmax * (Lmax + 2);

    % Partition Sg into sub-blocks:
    %   Sg = [ Gam  R
    %          T    S ]
    Gam = Sg(1:nPortmodes,              1:nPortmodes);              % Port-to-port
    R   = Sg(1:nPortmodes,              nPortmodes+1:end);          % Port-to-SWF
    T   = Sg(nPortmodes+1:end,          1:nPortmodes);              % SWF-to-port
    U   = (Sg(nPortmodes+1:end, nPortmodes+1:end) - eye(Nswf)) / 2; % Internal SWF block

    % Compute analytical embedding operator for the layered sphere
    omega = 2 * pi * freqs(i) * 1e9;   % Angular frequency [rad/s]
    [t, Phi, Psi, Rho] = utilities.GenSSO_ana( ... %#ok<NASGU>
        Lmax, ...
        r, ...
        [epsi; epsi], ...
        [mu;   mu  ], ...
        omega / c);

    % Compute modified port scattering matrix Γ~ using embedding relation
    %   Γ~ = Γ + 1/2 * R * Rho * (I - U * Rho)^(-1) * T
    Gam_new = Gam + R / 2 * Rho / (eye(Nswf) - U * Rho) * T;

    % Store the diagonal port reflection terms and the 4–5 coupling term
    Spara_new(i, :) = [ ...
        Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
        Gam_new(4,4), Gam_new(5,5), Gam_new(4,5) ...
    ];
end

%% Plot results: 20*log10(|Γ~|)
fig = figure; %#ok<NASGU>
set(gcf, 'Color', 'w');

axesHandle = axes; %#ok<LAXES>
hold(axesHandle, 'on');

% Convert to dB and plot all six traces
hPlot = plot(freqs, 20 * log10(abs(Spara_new)), 'LineWidth', 1);

% Axis properties and labels
set(axesHandle, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8, ...
    'Box', 'on', ...
    'LineWidth', 0.5, ...
    'GridLineStyle', ':', ...
    'YMinorTick', 'off', ...
    'MinorGridLineStyle', 'none', ...
    'GridColor', 'k');

xlim([freqs(1), freqs(end)]);
xticks(freqs);
xlabel('Frequency (GHz)', 'Interpreter', 'none');
ylabel('$20\log_{10}|\tilde{\mathbf{\Gamma}}|$', ...
    'Interpreter', 'latex');

grid(axesHandle, 'on');

legend(axesHandle, hPlot, ...
    {'$\tilde{\mathbf{\Gamma}}_{11}$', ...
     '$\tilde{\mathbf{\Gamma}}_{22}$', ...
     '$\tilde{\mathbf{\Gamma}}_{33}$', ...
     '$\tilde{\mathbf{\Gamma}}_{44}$', ...
     '$\tilde{\mathbf{\Gamma}}_{55}$', ...
     '$\tilde{\mathbf{\Gamma}}_{45}$'}, ...
    'NumColumns', 3, ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8, ...
    'Box', 'on', ...
    'Location', 'southwest');
