function T = buildDItable(DI)
% buildDItable
% Parse DI-field material definitions into a material-property table.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   T = buildDItable(DI)
%
%   INPUT:
%     DI : N×1 string array, e.g.
%           "DI: Dielectric1 : 0 : -1 :  :  : 5 :  :  :  : 0 : 1000"
%
%   OUTPUT:
%     T  : table with the columns:
%            - Label  : material name (string)
%            - EpsRel : relative permittivity ε_r
%            - Sigma  : conductivity (S/m)
%            - TanD   : loss tangent
%
%   RULES:
%     • Label:
%         If DI contains no name or the name is "0", the label is "".
%     • Slot #6  → relative permittivity (EpsRel), default = 1.
%     • Slot #8  → conductivity (Sigma), default = 0.
%     • Slot #10 → loss tangent (TanD), default = 0.
%
%   NOTES:
%     • Only dielectric materials are supported in this parsing routine.
%     • Missing or non-numeric values in DI are replaced by defaults.

    DI = string(DI(:));
    n  = numel(DI);

    Label  = strings(n,1);
    EpsRel = ones(n,1);    % default: vacuum permittivity
    Sigma  = zeros(n,1);   % default: no conductivity
    TanD   = zeros(n,1);   % default: no loss

    for i = 1:n
        line = strtrim(DI(i));
        if strlength(line) == 0
            continue;
        end

        % Split fields and strip the DI label
        parts = split(line, " : ");
        parts(1) = strtrim(erase(parts(1), "DI:"));

        % Ensure we have at least 11 fields
        if numel(parts) < 11
            parts(end+1:11) = "";
        end

        % --- Material label ---
        lab = strtrim(parts(1));
        if lab == "" || lab == "0"
            lab = "";
        end
        Label(i) = lab;

        % --- EpsRel (slot 6) ---
        v = str2double(strtrim(parts(6)));
        if ~isnan(v)
            EpsRel(i) = v;
        else
            EpsRel(i) = 1;
        end

        % --- Sigma (slot 8) ---
        v = str2double(strtrim(parts(8)));
        if ~isnan(v)
            Sigma(i) = v;
        else
            Sigma(i) = 0;
        end

        % --- TanD (slot 10) ---
        v = str2double(strtrim(parts(10)));
        if ~isnan(v)
            TanD(i) = v;
        else
            TanD(i) = 0;
        end
    end

    T = table(Label, EpsRel, Sigma, TanD);
end
