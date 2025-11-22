function portsInfo = annotatePortsFromAW(portTables, AW)
% annotatePortsFromAW
% Build a port-information structure array from AW lines.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   portsInfo = annotatePortsFromAW(portTables, AW)
%
%   This function parses AW definition lines and returns a structure array
%   describing each detected port. It does NOT modify or write back to
%   portTables.
%
%   INPUTS:
%     portTables : containers.Map whose keys are port names (char / string)
%                  and whose values are tables with at least the columns:
%                       - Label : string or char, point label
%                       - X, Y, Z: double, coordinates of the point
%     AW         : string/cellstr/char vector of lines containing "AW:"
%
%   OUTPUT:
%     portsInfo  : 1×N struct array with fields
%                  .portName    (char)     : port name
%                  .faceLabel   (string)   : 3rd field in AW (geometry face)
%                  .typeCode    (double)   : type code parsed from AW
%                  .typeName    (string)   : mapped type name (via code2name)
%                  .mediumLabel (string)   : material label (may be empty "")
%                  .ptLabs      (string[]) : point labels, in AW order
%                  .XYZ         (double)   : 3×N coordinate matrix (X;Y;Z)
%
%   NOTE:
%   - The mapping from type code to type name (code2name) is a placeholder
%     and should be adapted to your specific AW format.

    % Ensure AW is a column string array
    AW = string(AW(:));

    % Type-code to type-name mapping (adapt to your own definitions)
    code2name = containers.Map( ...
        {'1','2','3'}, ...
        {'Rectangular','Circular','Coaxial'} ...
    );

    % Pre-allocate output as empty struct array
    portsInfo = struct( ...
        'portName',    {}, ...
        'faceLabel',   {}, ...
        'typeCode',    {}, ...
        'typeName',    {}, ...
        'mediumLabel', {}, ...
        'ptLabs',      {}, ...
        'XYZ',         {} );

    % --- Main parsing loop over AW lines ---
    for i = 1:numel(AW)
        line = strtrim(AW(i));

        % Skip non-AW lines
        if ~startsWith(line, "AW:")
            continue;
        end

        % Split into "left part" and optional "Source" part
        seg  = split(line, "  ** ");
        left = strtrim(seg(1));

        % Fields after "AW:":
        %   toks(1) = "AW: <idx>"
        %   toks(2) = face label
        %   toks(3) = type code (string)
        %   toks(4) = medium label (optional, may be empty)
        %   toks(5..end) = remaining tokens / point labels
        toks = split(left, " : ");
        if numel(toks) < 3
            continue;
        end

        faceLabel   = strtrim(toks(2));
        typeCodeStr = strtrim(toks(3));
        typeCode    = str2double(typeCodeStr);

        if isKey(code2name, typeCodeStr)
            typeName = string(code2name(typeCodeStr));
        else
            typeName = "Unknown";
        end

        % Medium label (optional, can be empty string)
        if numel(toks) >= 4
            mediumLabel = string(strtrim(toks(4)));
            if mediumLabel == ""
                mediumLabel = "";
            end
        else
            mediumLabel = "";
        end

        % Collect point labels:
        % Prefer tokens from index 5 onward, which is more robust.
        ptLabs = strings(0, 1);
        for k = 5:numel(toks)
            tk = strtrim(toks(k));
            % Point labels are assumed to end with "_R<number>"
            if ~isempty(tk) && ~isempty(regexp(tk, '_R\d+$', 'once'))
                ptLabs(end+1, 1) = tk; %#ok<AGROW>
            end
        end

        % Backward compatibility:
        % If nothing found and token 4 looks like a point label, use it.
        if isempty(ptLabs) && numel(toks) >= 4
            tk = strtrim(toks(4));
            if ~isempty(tk) && ~isempty(regexp(tk, '_R\d+$', 'once'))
                ptLabs = tk;
            end
        end

        % If no point labels detected, skip this AW line
        if isempty(ptLabs)
            continue;
        end

        % Derive port name from first point label
        portName = char(localGetPortName(ptLabs(1)));

        % If the port name is not present in portTables, skip
        if ~isKey(portTables, portName)
            continue;
        end

        % Extract XYZ coordinates for all points (in AW order)
        T = portTables(portName); % Expected columns: Label, X, Y, Z

        % Find indices of AW point labels in the table
        [tfAW, idxInT] = ismember(ptLabs, T.Label);
        idxInT = idxInT(tfAW);

        if isempty(idxInT)
            XYZ          = zeros(3, 0);
            ptLabsFound  = strings(0, 1);
        else
            X = T.X(idxInT);
            Y = T.Y(idxInT);
            Z = T.Z(idxInT);
            XYZ         = [X.'; Y.'; Z.'];  % 3×N, columns correspond to points
            ptLabsFound = ptLabs(tfAW);
        end

        % Append to portsInfo
        portsInfo(end+1) = struct( ... %#ok<AGROW>
            'portName',    portName, ...
            'faceLabel',   string(faceLabel), ...
            'typeCode',    typeCode, ...
            'typeName',    string(typeName), ...
            'mediumLabel', string(mediumLabel), ...
            'ptLabs',      ptLabsFound, ...
            'XYZ',         XYZ );
    end
end

% --------- Local helper function ---------
function pn = localGetPortName(lab)
% localGetPortName
% Extract a port name from a point label, e.g.
%   "Port_C_Port_C_R3" -> "Port_C"
%   "PORT1_R2"         -> "PORT1"
%
% The extraction logic is:
%   1) Remove trailing "_R<number>"
%   2) If the remaining string is "<name>_<name>", return "<name>"
%   3) Otherwise, return the substring before the first underscore
%      (or the full string if no underscore is present).

    base = regexprep(string(lab), '_R\d+$', '');
    tok  = regexp(base, '^(.+)_\1$', 'tokens', 'once');
    if ~isempty(tok)
        pn = string(tok{1});
    else
        tok2 = regexp(base, '^([^_]+)', 'tokens', 'once');
        if ~isempty(tok2)
            pn = string(tok2{1});
        else
            pn = base;
        end
    end
end
