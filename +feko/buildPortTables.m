function portTables = buildPortTables(DP, varargin)
% buildPortTables  Parse DP strings into per-port point tables.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   portTables = buildPortTables(DP)
%   portTables = buildPortTables(DP, 'ZeroTol', value)
%
%   INPUTS:
%     DP       : N×1 string (or char/cellstr convertible to string), each line
%                describing a port point, e.g.
%                "DP: Port_C_Port_C_R3 : ... : X : Y : Z"
%     ZeroTol  : (optional, name-value) numerical threshold. Any |coord| <
%                ZeroTol is set to 0. Default: 1e-12.
%
%   OUTPUT:
%     portTables : containers.Map keyed by port name (char). Each value is a
%                  table with variables:
%                     - Label : original point label (string)
%                     - X, Y, Z : coordinates (double), sorted by R-index
%                                if present in the label.
%
%   The function:
%     * extracts a port name from each DP label (e.g. "Port_C" from
%       "Port_C_Port_C_R3"),
%     * groups rows by port name,
%     * sorts points within each port by their trailing "_R<number>" index
%       when available.

    % ----------------- Parse inputs -----------------
    p = inputParser;
    addParameter(p, 'ZeroTol', 1e-12, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    zeroTol = p.Results.ZeroTol;

    % Normalize DP to column string array
    DP = string(DP(:));
    n  = numel(DP);

    % ----------------- Preallocation -----------------
    labels = strings(n, 1);
    ports  = strings(n, 1);
    X      = nan(n, 1);
    Y      = nan(n, 1);
    Z      = nan(n, 1);
    Ridx   = nan(n, 1);

    % ----------------- Line parsing -----------------
    for i = 1:n
        parts = split(DP(i), " : ");
        if numel(parts) < 4
            % Not enough fields to contain coordinates; skip gracefully
            continue;
        end

        % "DP: <Label>" -> <Label>
        head = strtrim(parts(1));
        lab  = strtrim(erase(head, "DP:"));
        labels(i) = lab;

        % Coordinates assumed in the last three fields
        X(i) = str2double(strtrim(parts(end-2)));
        Y(i) = str2double(strtrim(parts(end-1)));
        Z(i) = str2double(strtrim(parts(end)));

        % --------- Extract port name ---------
        % 1) Remove trailing "_R<number>"
        base = regexprep(lab, '_R\d+$', '');
        % 2) Prefer form "<PortName>_<PortName>"
        tok = regexp(base, '^(.+)_\1$', 'tokens', 'once');
        if ~isempty(tok)
            ports(i) = string(tok{1});             % e.g. "Port_C"
        else
            % Fallback: use base; if that still looks too long, take prefix
            ports(i) = string(base);
            tok2 = regexp(lab, '^([^_]+)_', 'tokens', 'once');
            if ~isempty(tok2)
                ports(i) = string(tok2{1});
            end
        end

        % --------- Extract R index (if any) ---------
        rtok = regexp(lab, '_R(\d+)$', 'tokens', 'once');
        if ~isempty(rtok)
            Ridx(i) = str2double(rtok{1});
        end
    end

    % Zero-out very small coordinates
    X(abs(X) < zeroTol) = 0;
    Y(abs(Y) < zeroTol) = 0;
    Z(abs(Z) < zeroTol) = 0;

    % ----------------- Group by port -----------------
    up = unique(ports);
    portTables = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for k = 1:numel(up)
        pk   = up(k);
        mask = (ports == pk);

        lab_k = labels(mask);
        X_k   = X(mask);
        Y_k   = Y(mask);
        Z_k   = Z(mask);
        R_k   = Ridx(mask);

        % Sort by R index (if present). Points with NaN R-index are placed last.
        if any(~isnan(R_k))
            [~, ord] = sortrows([isnan(R_k), R_k]);  % numbered first, ascending R
        else
            ord = 1:numel(lab_k);
        end

        T = table(lab_k(ord), X_k(ord), Y_k(ord), Z_k(ord), ...
                  'VariableNames', {'Label', 'X', 'Y', 'Z'});

        portTables(char(pk)) = T;
    end
end
