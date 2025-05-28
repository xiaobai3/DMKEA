function [Offspring] = RealVariation_MutateRateByActive(Problem, ParentDec, OffMask)
% RealVariation_MutateRateByActive - Real-coded variation guided by binary mask
%
% This operator applies SBX crossover and polynomial mutation only on variables
% selected by a binary mask. Mutation parameters (disM, disC) are dynamically
% adjusted based on the progress ratio (delta).
%
% Developed by: Lidan Bai, 2024

    %% Pair parents
    Parent1 = ParentDec(1:floor(end/2), :);
    Parent2 = ParentDec(floor(end/2)+1 : floor(end/2)*2, :);

    %% Progress ratio
    delta = Problem.FE / Problem.maxFE;

    %% Parameter setting
    [N, D] = size(Parent1);
    [proC, ~, proM, ~] = deal(1, 20, 1, 20); 

    % Adaptive distribution indices
    disM_min = 20;
    disM_max = 100;
    disM = disM_min + (disM_max - disM_min) * delta^2;
    disC = disM;

    %% SBX Crossover
    beta = zeros(N, D);
    mu = rand(N, D);
    beta(mu <= 0.5) = (2 * mu(mu <= 0.5)).^(1 / (disC + 1));
    beta(mu > 0.5)  = (2 - 2 * mu(mu > 0.5)).^(-1 / (disC + 1));
    beta = beta .* (-1).^randi([0, 1], N, D);
    beta(rand(N, D) < 0.5) = 1;
    beta(repmat(rand(N, 1) > proC, 1, D)) = 1;

    Offspring = (Parent1 + Parent2) / 2 + beta .* (Parent1 - Parent2) / 2;
    Lower = repmat(Problem.lower, N, 1);
    Upper = repmat(Problem.upper, N, 1);

    %% Polynomial mutation (selective: only mutate active variables)
    d = sum(OffMask, 2);
    d = repmat(d, 1, D);
    Site = OffMask & rand(N, D) < proM ./ d;

    mu = rand(N, D);
    Offspring = min(max(Offspring, Lower), Upper);

    % Mutation where mu <= 0.5
    temp = Site & mu <= 0.5;
    Offspring(temp) = Offspring(temp) + (Upper(temp) - Lower(temp)) .* ...
        ((2 .* mu(temp) + (1 - 2 .* mu(temp)) .* ...
        (1 - (Offspring(temp) - Lower(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1);

    % Mutation where mu > 0.5
    temp = Site & mu > 0.5;
    Offspring(temp) = Offspring(temp) + (Upper(temp) - Lower(temp)) .* ...
        (1 - (2 .* (1 - mu(temp)) + 2 .* (mu(temp) - 0.5) .* ...
        (1 - (Upper(temp) - Offspring(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)));
end
