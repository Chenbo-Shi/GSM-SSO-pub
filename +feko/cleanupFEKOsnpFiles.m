function cleanupFEKOsnpFiles(preFEKOfileName)
% CLEANUPFEKOSNPFILES
%   Remove temporary FEKO S-parameter files generated during simulation.
%   This function deletes files matching:
%       <name>_SPara.s*p
%
% Notes
%   - Logic is unchanged; only documentation and copyright have been
%     updated to meet publication requirements.
%
% (c) 2025, Chenbo Shi, UESTC in China
% -------------------------------------------------------------------------

try
    delete([preFEKOfileName '_SPara.s*p']);
end

end
