function cleanupFEKOfarfieldFiles(preFEKOfileName)
% CLEANUPFEKOFARFIELDFILES
%   Remove temporary FEKO far-field files generated during simulation.
%   This function deletes all files matching:
%       <name>_FarField*.ffe
%       <name>_PortFarField*.ffe
%
% Notes
%   - Computational behavior unchanged; only documentation and copyright
%     have been updated for release.
%
% (c) 2025, Chenbo Shi, UESTC in China
% -------------------------------------------------------------------------

try
    delete([preFEKOfileName '_FarField*.ffe']);
end

try
    delete([preFEKOfileName '_PortFarField*.ffe']);
end

end
