function cleanupFEKOconfigfile()
% CLEANUPFEKOCONFIGFILE
%   Remove temporary FEKO runtime configuration files generated during
%   execution. This function deletes all files matching the pattern:
%       runfeko_tmp_configfile*.*
%
% Notes
%   - No logic has been altered; only documentation and copyright have been
%     updated for publication.
%
% (c) 2025, Chenbo Shi, UESTC in China
% -------------------------------------------------------------------------

try
    delete('runfeko_tmp_configfile*.*');
end

end
