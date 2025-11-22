function [F_th, F_ph, theta, phi, freq] = read_ffe_file(filename)
% READ_FFE_FILE  Read a FEKO .ffe far-field file into 3-D matrices.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [F_th, F_ph, theta, phi, freq] = READ_FFE_FILE(filename)
%
%   This function parses a FEKO .ffe far-field file and extracts only the
%   essential information: the column-name header and the '#Frequency' lines.
%   All other metadata is ignored.
%
%   OUTPUTS:
%     F_th(theta,phi,freq) - complex E_theta field samples
%     F_ph(theta,phi,freq) - complex E_phi   field samples
%     theta                - unique sorted polar angles (as in file)
%     phi                  - unique sorted azimuth angles (as in file)
%     freq                 - unique sorted frequencies (Hz)

    fid = fopen(filename, 'rt');
    if fid < 0
        error('read_ffe_file:FileOpenFailed', ...
              'Failed to open file: %s', filename);
    end
    c = onCleanup(@() fclose(fid)); 

    % Accumulators
    ThetaList = [];
    PhiList   = [];
    FreqList  = [];
    ReEthList = [];
    ImEthList = [];
    ReEphList = [];
    ImEphList = [];

    curFreq   = NaN;   % current frequency from "#Frequency" line
    colMap    = struct('theta',[], 'phi',[], ...
                       'reEth',[], 'imEth',[], 'reEph',[], 'imEph',[]);
    haveColMap = false;  % becomes true after a quoted-column header is parsed

    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        s = strtrim(tline);
        if isempty(s)
            continue;
        end

        % -----------------------------------------------------------------
        % 1) Header / comment lines: start with # / ## / ** (ignored except
        %    for #Frequency and quoted column-name line)
        % -----------------------------------------------------------------
        if startsWith(s, '#') || startsWith(s, '##') || startsWith(s, '**')
            % Parse frequency
            if contains(s, '#Frequency:')
                tok = regexp(s, '#Frequency:\s*([+\-]?\d+(\.\d+)?([Ee][+\-]?\d+)?)', ...
                             'tokens', 'once');
                if ~isempty(tok)
                    curFreq = str2double(tok{1});
                else
                    curFreq = NaN;
                end
            end

            % Parse column names if present (quoted fields)
            if contains(s, '"Theta"') || contains(s, '"Phi"') || ...
               contains(lower(s), 're(etheta)')
                names = regexp(s, '"([^"]+)"', 'tokens');
                names = cellfun(@(x) x{1}, names, 'UniformOutput', false);

                % Normalization for comparison
                normName = @(x) lower(regexprep(x, '\s+', ''));

                colMap = struct('theta',[], 'phi',[], ...
                                'reEth',[], 'imEth',[], 'reEph',[], 'imEph',[]);
                for k = 1:numel(names)
                    nm = normName(names{k});
                    if strcmp(nm, 'theta')
                        colMap.theta = k;
                    elseif strcmp(nm, 'phi')
                        colMap.phi = k;
                    elseif any(strcmp(nm, {'re(etheta)', 're-etheta', 're_etheta'}))
                        colMap.reEth = k;
                    elseif any(strcmp(nm, {'im(etheta)', 'im-etheta', 'im_etheta'}))
                        colMap.imEth = k;
                    elseif any(strcmp(nm, {'re(ephi)', 're-ephi', 're_ephi'}))
                        colMap.reEph = k;
                    elseif any(strcmp(nm, {'im(ephi)', 'im-ephi', 'im_ephi'}))
                        colMap.imEph = k;
                    end
                end
                haveColMap = all(structfun(@(v) ~isempty(v), colMap));
            end
            continue;
        end

        % -----------------------------------------------------------------
        % 2) Data lines (numeric)
        % -----------------------------------------------------------------
        nums = sscanf(s, '%f');
        if isempty(nums)
            continue;
        end

        % Fallback: if no explicit column map, assume FEKO standard order:
        % [Theta, Phi, Re(Etheta), Im(Etheta), Re(Ephi), Im(Ephi), ...]
        if ~haveColMap
            if numel(nums) < 6
                continue;
            end
            th  = nums(1);
            ph  = nums(2);
            reT = nums(3);
            imT = nums(4);
            reP = nums(5);
            imP = nums(6);
        else
            need = [colMap.theta, colMap.phi, ...
                    colMap.reEth, colMap.imEth, ...
                    colMap.reEph, colMap.imEph];
            if numel(nums) < max(need)
                continue;
            end
            th  = nums(colMap.theta);
            ph  = nums(colMap.phi);
            reT = nums(colMap.reEth);
            imT = nums(colMap.imEth);
            reP = nums(colMap.reEph);
            imP = nums(colMap.imEph);
        end

        ThetaList(end+1,1) = th;   %#ok<AGROW>
        PhiList(end+1,1)   = ph;   %#ok<AGROW>
        FreqList(end+1,1)  = curFreq; %#ok<AGROW>

        ReEthList(end+1,1) = reT;  %#ok<AGROW>
        ImEthList(end+1,1) = imT;  %#ok<AGROW>
        ReEphList(end+1,1) = reP;  %#ok<AGROW>
        ImEphList(end+1,1) = imP;  %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % 3) Assemble 3-D matrices (theta, phi, freq)
    % ---------------------------------------------------------------------
    if isempty(FreqList)
        error('read_ffe_file:NoDataParsed', ...
              'No far-field data found in file: %s', filename);
    end

    theta = unique(ThetaList, 'sorted');
    phi   = unique(PhiList,   'sorted');
    freq  = unique(FreqList,  'sorted');

    nTh = numel(theta);
    nPh = numel(phi);
    nF  = numel(freq);

    F_th = complex(nan(nTh, nPh, nF), nan(nTh, nPh, nF));
    F_ph = complex(nan(nTh, nPh, nF), nan(nTh, nPh, nF));

    [~, thIdx] = ismember(ThetaList, theta);
    [~, phIdx] = ismember(PhiList,   phi);
    [~, frIdx] = ismember(FreqList,  freq);

    Eth = ReEthList + 1i * ImEthList;
    Eph = ReEphList + 1i * ImEphList;

    lin = sub2ind([nTh, nPh, nF], thIdx, phIdx, frIdx);
    F_th(lin) = Eth;
    F_ph(lin) = Eph;
end
