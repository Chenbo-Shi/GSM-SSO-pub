function [normLegendreSing, normLegendreSingDer, normLegendre] = ...
    legendreComponents(degreeL, orderM, cosTheta)
% legendreComponents
% Computes normalized associated Legendre polynomials, their singularity-
% treated versions and derivatives for given degree L and order M.
%
%   [normLegendreSing, normLegendreSingDer, normLegendre] = ...
%                       legendreComponents(degreeL, orderM, cosTheta)
%
%   INPUTS
%     degreeL : degree l of the Legendre polynomial, double [N x 1]
%     orderM  : order m corresponding to each degreeL entry, double [N x 1]
%     cosTheta: cosine of the polar angle theta at integration points,
%               double [1 x K] or [K x 1]
%
%   OUTPUTS
%     normLegendreSing    : (m / sin(theta)) * P_l^m(cosTheta) with
%                           singularity treatment, double [N x K]
%     normLegendreSingDer : derivative w.r.t. theta of the singularity
%                           treated polynomial, double [N x K]
%     normLegendre        : normalized associated Legendre polynomial
%                           P_l^m(cosTheta), double [N x K]
%
% (c) 2025, Chenbo Shi, UESTC in China

%% Input shaping
degreeL = degreeL(:);
orderM  = orderM(:);
nModes  = numel(degreeL);
nTheta  = numel(cosTheta);

%% Pre-calculation of constants (vectorized over modes)
singConst1 = 0.5 * sqrt((2 .* degreeL + 1) ./ (2 .* degreeL + 3));
singConst2 = sqrt((degreeL + orderM + 2) .* (degreeL + orderM + 1));
singConst3 = sqrt((degreeL - orderM + 2) .* (degreeL - orderM + 1));

singConstDer1 = sqrt((degreeL + 1) .* degreeL);
singConstDer2 = sqrt((degreeL + orderM) .* (degreeL - orderM + 1));

maxDegreeL = max(degreeL);

%% Memory allocation
% normalized Legendre with singularity solution incorporated
normLegendreSing          = nan(nModes, nTheta);
normLegendreSingMplus1    = normLegendreSing;
normLegendreSingMminus1   = normLegendreSing;

% derivative of normalized Legendre with singularity solution incorporated
normLegendreSingMminus1_D = normLegendreSing;
normLegendreSingM1        = normLegendreSing;

% normalized Legendre
normLegendre              = normLegendreSing;

%% Derivative of normalized Legendre via recursive relations
temp = legendre(1, cosTheta, 'norm');  % l = 1

for thisDegreeL = 1:maxDegreeL
    % Normalized Legendre for l = thisDegreeL + 1
    temp1 = legendre(thisDegreeL + 1, cosTheta, 'norm');
    
    idxL = (degreeL == thisDegreeL);
    if ~any(idxL)
        % No entries with this degree; shift temp and continue
        temp = temp1;
        continue;
    end
    
    % Orders corresponding to this degree
    mThis = orderM(idxL);
    
    % --- P_{l+1}^{m+1} ---
    Mplus1 = temp1(mThis + 2, :);          % row index m+1 -> m+2 (1-based)
    
    % --- P_{l+1}^{m-1} ---
    Mminus1 = zeros(size(Mplus1));
    tempIndex = mThis;
    tempIndex(tempIndex == 0) = [];        % m = 0 has no m-1 term
    if ~isempty(tempIndex)
        Mminus1(mThis ~= 0, :) = temp1(tempIndex, :);
    end
    
    normLegendreSingMplus1(idxL, :)  = Mplus1;
    normLegendreSingMminus1(idxL, :) = Mminus1;
    
    % --- derivative-related part using l = thisDegreeL ---
    Mminus1 = zeros(sum(idxL), nTheta);
    if ~isempty(tempIndex)
        Mminus1(mThis ~= 0, :) = temp(tempIndex, :);        
    end
    
    normLegendreSingMminus1_D(idxL, :) = Mminus1;
    
    % Special handling for m = 0
    isM0 = (orderM == 0) & (degreeL == thisDegreeL);
    if any(isM0)
        % P_l^{1} for derivative relation, row index 2
        normLegendreSingM1(isM0, :) = repmat(temp(2, :), sum(isM0), 1);
    end
    
    % Store normalized P_l^m for this degree
    normLegendre(idxL, :) = temp(orderM(idxL) + 1, :);
    
    % Shift for next iteration: temp is now l = thisDegreeL + 1
    temp = temp1;
end

%% Assemble singularity-treated components
normLegendreSing = singConst1 .* ...
    (singConst2 .* normLegendreSingMplus1 + ...
     singConst3 .* normLegendreSingMminus1);

normLegendreSing(orderM == 0, :) = 0;

%% Derivative of singularity-treated Legendre
% cosTheta.' is [nTheta x 1], implicit expansion to [N x nTheta]
normLegendreSingDer = singConstDer2 .* normLegendreSingMminus1_D - ...
    cosTheta.' .* normLegendreSing;

% Special case for m = 0
normLegendreSingDer(orderM == 0, :) = ...
    -singConstDer1(orderM == 0) .* normLegendreSingM1(orderM == 0, :);
end
