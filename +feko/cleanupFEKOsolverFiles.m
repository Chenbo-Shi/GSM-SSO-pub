function cleanupFEKOsolverFiles(preFEKOfileName, skipTheseFiles)
% CLEANUPFEKOSOLVERFILES
%   Remove temporary FEKO solver files created during preprocessing and
%   solution stages of a FEKO run.
%
% Inputs
%   preFEKOfileName ~ base name of the FEKO files to delete
%   skipTheseFiles  ~ cell array of file extensions to skip (e.g. {'.pre'})
%
% Notes
%   - Only documentation and copyright have been updated.
%   - Computational behavior remains identical to the original.
%
% (c) 2025, Chenbo Shi, UESTC in China
% -------------------------------------------------------------------------

nInputs = nargin;
if nInputs < 2
    skipTheseFiles = {};
end

cleanup = {'.bof', '.fek', '.out', '.str', '.pre'};

for thisFile = find(~ismember(cleanup, skipTheseFiles))
    try
        delete([preFEKOfileName cleanup{thisFile}]);
    end
end

end
