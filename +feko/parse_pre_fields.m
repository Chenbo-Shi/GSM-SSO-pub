function [SF, IN, DP, EG, PS, CG, AW, DI, FP] = parse_pre_fields(prePath)
% PARSE_PRE_FIELDS  Extract selected FEKO cards from a .pre file.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [SF, IN, DP, EG, PS, CG, AW, DI, FP] = PARSE_PRE_FIELDS(prePath)
%
%   Inputs
%   -------
%   prePath : char or string
%       Path to the .pre file.
%
%   Outputs
%   -------
%   SF : char
%       First line starting with "SF:" (empty if not found).
%   IN : string column vector
%       All lines starting with "IN" (may contain multiple entries).
%   DP : string column vector
%       All lines starting with "DP:".
%   EG : char
%       First line starting with "EG:" (empty if not found).
%   PS : char
%       First line starting with "PS:" (empty if not found).
%   CG : char
%       First line starting with "CG:" (empty if not found).
%   AW : string column vector
%       All lines starting with "AW:" (may contain multiple entries).
%   DI : string column vector
%       All lines starting with "DI:" (may contain multiple entries).
%   FP : string column vector
%       All lines starting with "FP:" (may contain multiple entries).
%
%   Notes
%   -----
%   - Comment lines starting with "**" are ignored.
%   - "IN" cards may have leading spaces and no colon.
%   - Only the first occurrence of SF, EG, PS, CG is returned.

    % Read file content and split into lines
    txt   = fileread(prePath);
    lines = regexp(txt, '\r\n|\n|\r', 'split');

    % Initialize outputs
    SF = '';
    EG = '';
    PS = '';
    CG = '';

    IN = strings(0, 1);
    DP = strings(0, 1);
    AW = strings(0, 1);
    DI = strings(0, 1);
    FP = strings(0, 1);

    % Parse line by line
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if isempty(ln) || startsWith(ln, '**')
            continue;  % skip empty or comment lines
        end

        % --- SF card (first occurrence only) ---
        if startsWith(ln, 'SF:', 'IgnoreCase', true)
            if isempty(SF)
                SF = ln;
            end
            continue;
        end

        % --- IN cards (may have leading spaces, no colon) ---
        if startsWith(regexprep(ln, '^\s+', ''), 'IN', 'IgnoreCase', true)
            IN(end+1, 1) = string(ln); %#ok<AGROW>
            continue;
        end

        % --- DP cards ---
        if startsWith(ln, 'DP:', 'IgnoreCase', true)
            DP(end+1, 1) = string(ln); %#ok<AGROW>
            continue;
        end

        % --- EG, PS, CG: first occurrence only ---
        if startsWith(ln, 'EG:', 'IgnoreCase', true)
            if isempty(EG)
                EG = ln;
            end
            continue;
        end

        if startsWith(ln, 'PS:', 'IgnoreCase', true)
            if isempty(PS)
                PS = ln;
            end
            continue;
        end

        if startsWith(ln, 'CG:', 'IgnoreCase', true)
            if isempty(CG)
                CG = ln;
            end
            continue;
        end

        % --- AW, DI, FP cards (all occurrences) ---
        if startsWith(ln, 'AW:', 'IgnoreCase', true)
            AW(end+1, 1) = string(ln); %#ok<AGROW>
            continue;
        end

        if startsWith(ln, 'DI:', 'IgnoreCase', true)
            DI(end+1, 1) = string(ln); %#ok<AGROW>
            continue;
        end

        if startsWith(ln, 'FP:', 'IgnoreCase', true)
            FP(end+1, 1) = string(ln); %#ok<AGROW>
            continue;
        end
    end

    % Convert scalar string outputs to char for backward compatibility
    SF = char(SF);
    EG = char(EG);
    PS = char(PS);
    CG = char(CG);
end
