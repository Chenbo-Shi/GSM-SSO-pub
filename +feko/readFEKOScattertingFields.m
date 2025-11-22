function Scolumn = readFEKOScattertingFields(preFEKOfileName, nFarFields, thisPW)
% readFEKOScattertingFields  Read FEKO .ffe far-field files and assemble one
% scattering-column vector (for all far-field directions).
% (c) 2025, Chenbo Shi, UESTC in China
%
%   Scolumn = readFEKOScattertingFields(preFEKOfileName, nFarFields, thisPW)
%
%   Inputs:
%     preFEKOfileName : base name of the FEKO project (without suffix)
%     nFarFields      : number of far-field directions (number of .ffe files)
%     thisPW          : optional plane-wave identifier suffix (e.g. '_PW1'),
%                       default is '' (empty)
%
%   Output:
%     Scolumn         : 2*nFarFields × Nfreq complex column block
%                       [Ftheta; Fphi], where each row corresponds to one
%                       far-field direction and columns correspond to
%                       frequencies.

    if nargin < 2
        error('readFEKOScattertingFields:NotEnoughInputs', ...
            'preFEKOfileName and nFarFields must be specified.');
    end
    if nargin < 3 || isempty(thisPW)
        thisPW = '';
    end

    % Base name of FEKO far-field files
    farfieldFileBase = [preFEKOfileName '_FarField'];

    Fth = [];   % will be allocated after first file is read
    Fph = [];
    nFreq = [];

    % Loop over all far-field directions
    for n = 1:nFarFields
        thisFileName = [farfieldFileBase num2str(n) thisPW '.ffe'];

        % read_ffe_file returns 3-D arrays: (theta, phi, freq)
        [F_th_this, F_ph_this, ~, ~, ~] = feko.read_ffe_file(thisFileName);

        % Use only the first (theta, phi) sample; keep all frequencies
        thisFth = squeeze(F_th_this(1, 1, :)).';
        thisFph = squeeze(F_ph_this(1, 1, :)).';

        if isempty(nFreq)
            % First file: determine number of frequencies and allocate arrays
            nFreq = numel(thisFth);
            Fth   = zeros(nFarFields, nFreq);
            Fph   = zeros(nFarFields, nFreq);
        else
            % Consistency check
            if numel(thisFth) ~= nFreq
                error('readFEKOScattertingFields:FreqCountMismatch', ...
                    'File %s has a different number of frequencies.', thisFileName);
            end
        end

        Fth(n, :) = thisFth;
        Fph(n, :) = thisFph;
    end

    % Arrange data as a single column block (scattering dyadic column)
    % First nFarFields rows: Ftheta, next nFarFields rows: Fphi
    Scolumn = [Fth; Fph];
end
