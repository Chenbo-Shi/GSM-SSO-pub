function [scale, unitStr] = getScaleFromSF(SF)
% getScaleFromSF
%   Determine geometric scale factor and unit string from an SF line.
% (c) 2025, Chenbo Shi, UESTC in China
%
%   [scale, unitStr] = getScaleFromSF(SF)
%
%   INPUT:
%     SF : string or cellstr, e.g. "SF: 1 : : : : : 0.001"
%
%   OUTPUT:
%     scale   : numeric scale factor converting model units to meters
%     unitStr : detected unit string (e.g. "m", "mm", "cm", ...)

    if isempty(SF)
        scale   = 1.0;
        unitStr = "m";
        return
    end

    SF = string(SF);
    parts = split(SF, " : ");
    if numel(parts) < 2
        scale   = 1.0;
        unitStr = "m";
        return
    end

    % last field is assumed to be the scale value
    val = str2double(strtrim(parts(end)));
    if isnan(val) || val <= 0
        val = 1.0;
    end

    scale = val;

    % match common units
    switch round(val, 6)
        case 1
            unitStr = "m";
        case 1e-3
            unitStr = "mm";
        case 1e-2
            unitStr = "cm";
        case 0.0254
            unitStr = "in";
        case 0.3048
            unitStr = "ft";
        case 1e-6
            unitStr = "µm";
        otherwise
            unitStr = sprintf("scale=%.6g m/unit", val);
    end
end
