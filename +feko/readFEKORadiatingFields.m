function Rcolumn = readFEKORadiatingFields(preFEKOfileName, nFarFields, thismode)
% READFEKORADIATINGFIELDS  Read FEKO port far-field .ffe files and assemble one column.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   Rcolumn = readFEKORadiatingFields(preFEKOfileName, nFarFields, thismode)
%
%   This function reads a sequence of FEKO far-field files corresponding to
%   different port excitations / far-field directions and stacks their
%   theta- and phi-polarized far fields into a single column vector.
%
%   INPUTS:
%     preFEKOfileName : base name of the FEKO model (without suffix)
%     nFarFields      : number of far-field files to read
%     thismode        : optional string appended to the file name
%                       (e.g., '_1', '_2'; default is '')
%
%   FILE NAME PATTERN:
%     <preFEKOfileName>_PortFarField<n><thismode>.ffe
%     where n = 1, 2, ..., nFarFields
%
%   OUTPUT:
%     Rcolumn : stacked far-field samples as one column vector:
%               [Fth(:); Fph(:)], where
%               Fth and Fph are size [nFarFields × nFreq].

    if nargin < 3
        thismode = '';
    end

    baseName = [preFEKOfileName '_PortFarField'];

    % Read first file to determine frequency dimension and preallocate
    firstFile = [baseName '1' thismode '.ffe'];
    [F_th_first, F_ph_first, ~, ~, ~] = feko.read_ffe_file(firstFile);
    % We only use the sample at (theta(1), phi(1), :)
    Fth = zeros(nFarFields, size(F_th_first, 3));
    Fph = zeros(nFarFields, size(F_ph_first, 3));
    Fth(1, :) = reshape(F_th_first(1, 1, :), 1, []);
    Fph(1, :) = reshape(F_ph_first(1, 1, :), 1, []);

    % Remaining files
    for n = 2:nFarFields
        thisFileName = [baseName num2str(n) thismode '.ffe'];
        [F_th_this, F_ph_this, ~, ~, ~] = feko.read_ffe_file(thisFileName);

        % Basic size check (frequency dimension must match)
        if size(F_th_this, 3) ~= size(Fth, 2)
            error('Frequency dimension mismatch in file: %s', thisFileName);
        end

        Fth(n, :) = reshape(F_th_this(1, 1, :), 1, []);
        Fph(n, :) = reshape(F_ph_this(1, 1, :), 1, []);
    end

    % Arrange data as a single column vector: [Fth(:); Fph(:)]
    Rcolumn = [Fth(:); Fph(:)];
end
