function [Offspring] = RealVariation_MutateRateFixed(Problem, ParentDec, OffMask)
% RealVariation_MutateRateFixed - SBX + Polynomial mutation for active variables
%
% This operator applies real-coded crossover and polynomial mutation to the
% active variables (OffMask==1). Mutation rate is fixed to 1/D.
% The variation is scaled by the current optimization progress (delta).
%
% Written by Lidan Bai, 2024

    %% Split parents
    Parent1 = ParentDec(1:floor(end/2), :);
    Parent2 = ParentDec(floor(end/2)+1:floor(end/2)*2, :);
    delta   = Problem.FE / Problem.maxFE;  % Normalized progress

    %% Parameter setting
    [N, D]           = size(Parent1);
    [proC, ~, proM, ~] = deal(1, 20, 1, 20); 
    disM_min = 20;
    disM_max = 100;
    disM     = disM_min + (disM_max - disM_min) * delta^2;
    disC     = disM;

    %% Simulated Binary Crossover (SBX)
    beta = zeros(N, D);
    mu   = rand(N, D);
    beta(mu <= 0.5) = (2 * mu(mu <= 0.5)).^(1 / (disC + 1));
    beta(mu > 0.5)  = (2 - 2 * mu(mu > 0.5)).^(-1 / (disC + 1));
    beta = beta .* (-1).^randi([0, 1], N, D);
    beta(rand(N, D) < 0.5) = 1;
    beta(repmat(rand(N, 1) > proC, 1, D)) = 1;
    Offspring = (Parent1 + Parent2) / 2 + beta .* (Parent1 - Parent2) / 2;

    %% Bound handling
    Lower = repmat(Problem.lower, N, 1);
    Upper = repmat(Problem.upper, N, 1);
    Offspring = min(max(Offspring, Lower), Upper);

    %% Polynomial Mutation (fixed rate for each variable: 1/D)
    Site = OffMask & rand(N, D) < proM / D;
    mu   = rand(N, D);

    % Lower half mutation
    temp = Site & mu <= 0.5;
    Offspring(temp) = Offspring(temp) + ...
        (Upper(temp) - Lower(temp)) .* ...
        ((2 .* mu(temp) + (1 - 2 .* mu(temp)) .* ...
        (1 - (Offspring(temp) - Lower(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1);

    % Upper half mutation
    temp = Site & mu > 0.5;
    Offspring(temp) = Offspring(temp) + ...
        (Upper(temp) - Lower(temp)) .* ...
        (1 - (2 .* (1 - mu(temp)) + 2 .* (mu(temp) - 0.5) .* ...
        (1 - (Upper(temp) - Offspring(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)));
end
