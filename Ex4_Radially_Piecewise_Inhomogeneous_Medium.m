%% Radially Piecewise Inhomogeneous Medium (Sec. V-D)
% (c) 2025, Chenbo Shi, UESTC in China
%
% This script compares:
%   1) An ODE-based evaluation of the spherical shell operator (radially
%      inhomogeneous, continuous medium).
%   2) A multi-layer piecewise-homogeneous approximation of the same medium.
%
% The generalized scattering matrix Sg of a 5-mode horn antenna in free space
% must be precomputed first by running the FEKO–MATLAB co-simulation script
% (Ex0). The resulting data are expected in:
%   preserved_data/CoSim_Ex0_Horn_Liang_Script/Horn_Liang_Script_data(i).mat

%% Problem definition
freqs = 3.2 : 0.05 : 3.8;     % Frequency range [GHz]
nfreq = numel(freqs);
c     = 299792458;            % Speed of light [m/s]

% Radii of concentric spherical layers [m] (from inner to outer).
% The last "inf" represents free space.
r = [150 165 180 inf] * 1e-3;

% Radially inhomogeneous, anisotropic material parameters (continuous in r)
epsi_t = @(z)[1 5*tan(pi/5./z)      2 + log(2./z - 5) 1];
epsi_r = @(z)[1 1 + exp(2*sin(4./z))         1./z     1];
mu_t   = @(z)[1 1 1 1];
mu_r   = @(z)[1 1 1 1];

% Location of precomputed free-space generalized scattering matrices Sg
dataDir  = fullfile('preserved_data', 'CoSim_Ex0_Horn_Liang_Script');
baseName = 'Horn_Liang_Script_data';

%% ODE-based SSO operator for radially inhomogeneous medium
% Spara_new: selected entries of the modified reflection matrix Γ~
% Columns: [Γ11 Γ22 Γ33 Γ44 Γ55 Γ45]
Spara_new = zeros(nfreq, 6);

for i = 1:nfreq
    f = freqs(i);  % Frequency in GHz
    fprintf('ODE-SSO: processing f = %.3f GHz\n', f);

    % Load precomputed generalized scattering matrix Sg in free space
    dataFile = fullfile(dataDir, sprintf('%s(%d).mat', baseName, i));
    load(dataFile, 'Sg', 'nPortmodes', 'Lmax'); %#ok<LOAD>

    % Spherical-wave sub-block size (consistent with Lmax)
    Nfwav = 2 * Lmax * (Lmax + 2);
    Nmods = nPortmodes;

    % Partition Sg into sub-blocks:
    %   Sg = [ Gam  R
    %          T    S ]
    Gam = Sg(1:Nmods,             1:Nmods);              % Port-to-port
    R   = Sg(1:Nmods,             Nmods+1:Nmods+Nfwav);  % Port-to-SWF
    T   = Sg(Nmods+1:Nmods+Nfwav, 1:Nmods);              % SWF-to-port
    U   = (Sg(Nmods+1:Nmods+Nfwav, Nmods+1:Nmods+Nfwav) - eye(Nfwav)) / 2;

    % ODE-based spherical shell operator for radially inhomogeneous medium
    k0 = 2*pi*f*1e9 / c;
    [t, Phi, Psi, Rho] = utilities.GenSSO_ode_parallel( ... %#ok<ASGLU>
        Lmax, r, epsi_t, epsi_r, mu_t, mu_r, k0);

    % Modified reflection matrix in the presence of the medium
    Gam_new = Gam + R/2 * Rho / (eye(Nfwav) - U * Rho) * T;

    % Store selected entries
    Spara_new(i, :) = [Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
                       Gam_new(4,4), Gam_new(5,5), Gam_new(4,5)];
end

%% Multi-layer piecewise-homogeneous approximation

N = 21;

% Number of original spherical interfaces
Nlay = numel(r);

% ----- 150 mm – 165 mm: discretise into N-1 homogeneous shells -----
rL1 = linspace(r(1), r(2), N);
rL1(1) = [];                 % remove inner interface, keep sublayer radii
nL1   = numel(rL1);

et1 = zeros(1, nL1);
er1 = zeros(1, nL1);
mt1 = zeros(1, nL1);
mr1 = zeros(1, nL1);

for idx = 1:nL1
    rho   = rL1(idx);
    etVal = epsi_t(rho);     % returns a row vector [et_layer1 et_layer2 et_layer3 et_layer4]
    erVal = epsi_r(rho);
    mtVal = mu_t(rho);
    mrVal = mu_r(rho);

    % Take the parameters of the 2nd shell in this radial interval
    et1(idx) = etVal(2);
    er1(idx) = erVal(2);
    mt1(idx) = mtVal(2);
    mr1(idx) = mrVal(2);
end

% ----- 165 mm – 180 mm: discretise into N-1 homogeneous shells -----
rL2 = linspace(r(2), r(3), N);
rL2(1) = [];
nL2   = numel(rL2);

et2 = zeros(1, nL2);
er2 = zeros(1, nL2);
mt2 = zeros(1, nL2);
mr2 = zeros(1, nL2);

for idx = 1:nL2
    rho   = rL2(idx);
    etVal = epsi_t(rho);     % row vector as above
    erVal = epsi_r(rho);
    mtVal = mu_t(rho);
    mrVal = mu_r(rho);

    % Take the parameters of the 3rd shell in this radial interval
    et2(idx) = etVal(3);
    er2(idx) = erVal(3);
    mt2(idx) = mtVal(3);
    mr2(idx) = mrVal(3);
end

% ----- Assemble piecewise-homogeneous spherical stack -----
r  = [r(1), rL1, rL2, inf];  % new interfaces (inner radius, sublayers, outer radius, free space)
et = [1,    et1, et2, 1];    % transverse permittivity for each region
er = [1,    er1, er2, 1];    % radial permittivity
mt = [1,    mt1, mt2, 1];    % transverse permeability
mr = [1,    mr1, mr2, 1];    % radial permeability

Nlay = numel(r);             % updated number of layers

% Spara_multi: selected entries of Γ~ using the multi-layer approximation
Spara_multi = zeros(nfreq, 6);

for i = 1:nfreq
    f = freqs(i);
    fprintf('Multi-layer SSO: processing f = %.3f GHz\n', f);

    % Load precomputed generalized scattering matrix Sg in free space
    dataFile = fullfile(dataDir, sprintf('%s(%d).mat', baseName, i));
    load(dataFile, 'Sg', 'nPortmodes', 'Lmax'); %#ok<LOAD>

    % Spherical-wave sub-block size
    Nfwav = 2 * Lmax * (Lmax + 2);
    Nmods = nPortmodes;

    % Partition Sg into sub-blocks:
    Gam = Sg(1:Nmods,             1:Nmods);
    R   = Sg(1:Nmods,             Nmods+1:Nmods+Nfwav);
    T   = Sg(Nmods+1:Nmods+Nfwav, 1:Nmods);
    U   = (Sg(Nmods+1:Nmods+Nfwav, Nmods+1:Nmods+Nfwav) - eye(Nfwav)) / 2;

    % Analytical spherical shell operator for the multi-layer approximation
    k0 = 2*pi*f*1e9 / c;
    [t, Phi, Psi, Rho] = utilities.GenSSO_ana( ... %#ok<ASGLU>
        Lmax, r, [et; er], [mt; mr], k0);

    % Modified reflection matrix
    Gam_new = Gam + R/2 * Rho / (eye(Nfwav) - U * Rho) * T;

    Spara_multi(i, :) = [Gam_new(1,1), Gam_new(2,2), Gam_new(3,3), ...
                         Gam_new(4,4), Gam_new(5,5), Gam_new(4,5)];
end

%% Plot comparison of ODE vs. multi-layer approximation
fig = figure;

hold on;
% Overall ODE curves (black, all entries stacked)
h0 = plot(freqs, 20*log10(abs(Spara_new)), 'k', 'LineWidth', 1);
% Per-entry colored curves (ODE result)
h1 = plot(freqs, 20*log10(abs(Spara_new)), 'LineWidth', 1);
% Multi-layer approximation (markers)
h2 = plot(freqs, 20*log10(abs(Spara_multi)), 'o', ...
    'Color', '#636363', 'MarkerSize', 2, ...
    'MarkerIndices', 1:numel(freqs));

set(gca, 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'LineWidth', 0.5, ...
    'GridLineStyle', ':', 'yMinorTick', 'off', ...
    'MinorGridLineStyle', 'none', 'GridColor', 'k');

xlim([freqs(1), freqs(end)]);
xticks(freqs);
ylabel('$20\log_{10}|\tilde{\mathbf{\Gamma}}|$', 'Interpreter', 'latex');
xlabel('Frequency (GHz)');
grid on;

% First legend: individual matrix elements
ax1 = gca;
ax2 = axes('Parent', fig, 'Position', ax1.Position, 'Visible', 'off');

leg1 = legend(ax1, h1, { ...
    '$\tilde{\mathbf{\Gamma}}_{11}$', ...
    '$\tilde{\mathbf{\Gamma}}_{22}$', ...
    '$\tilde{\mathbf{\Gamma}}_{33}$', ...
    '$\tilde{\mathbf{\Gamma}}_{44}$', ...
    '$\tilde{\mathbf{\Gamma}}_{55}$', ...
    '$\tilde{\mathbf{\Gamma}}_{45}$'}, ...
    'NumColumns', 3, 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'Location', 'southwest'); 

% Second legend: method labels (ODE vs. multi-layer)
leg2 = legend(ax2, [h0(1), h2(1)], {'ODE', '10-Layers'}, ...
    'FontSize', 8, 'FontName', 'Times New Roman', ...
    'NumColumns', 1, 'Box', 'on', 'Location', 'southeast'); %#ok<NASGU>

set(gcf, 'Color', 'w');
set(ax1, 'GridLineStyle', ':', 'GridColor', 'k');

%% Local helper: extract i-th component from a vector-valued function handle
function f = get(fhandle, i)
% GET  Wrap a vector-valued function handle to extract its i-th component.
%   f = GET(fhandle, i) returns a new handle f(z) that evaluates the
%   original handle fhandle(z) and returns only its i-th component.
    f = @(z) subsref(fhandle(z), substruct('()', {{i}}));
end
