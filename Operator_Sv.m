function [OffDec, OffMask] = Operator_Sv(Problem, ParentDec, ParentMask, sv)
% Operator_Sv - Binary variation guided by statistical vector (sv)
%
% This operator reuses the sv-guided binary crossover and mutation design from MSKEA.
% The real-valued variation part is customized: only active variables are mutated
% using a fixed mutation rate (1/D), implemented in RealVariation_MutateRateFixed.
%
% Binary part: reproduced from MSKEA
% Real part:  implemented by Lidan Bai, 2024

    %% Parameter setting
    [N, D] = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2, :);
    Parent2Mask = ParentMask(N/2+1:end, :);

    %% Binary crossover guided by sv
    OffMask = Parent1Mask;
    rate0 = sv;          % Probability that 0 → 1
    rate1 = 1 - rate0;   % Probability that 1 → 0

    for i = 1:N/2
        diff = find(Parent1Mask(i,:) ~= Parent2Mask(i,:));
        temp_rate1 = rate1(diff);
        temp_rate0 = rate0(diff);
        rate = zeros(1, length(diff));
        rate(logical(OffMask(i, diff)))  = temp_rate1(logical(OffMask(i, diff)));
        rate(logical(~OffMask(i, diff))) = temp_rate0(logical(~OffMask(i, diff)));
        exchange = rand(1, length(diff)) < rate;
        OffMask(i, diff(exchange)) = ~OffMask(i, diff(exchange));
    end

    %% Binary mutation guided by sv
    Mutation_p = 1 / D;
    Mu_exchange = rand(N/2, D) < Mutation_p;
    for i = 1:N/2
        if any(Mu_exchange(i, :))
            subscript = find(Mu_exchange(i, :) == 1);
            rate = zeros(1, numel(subscript));
            rate(logical(OffMask(i, subscript)))  = rate1(subscript(logical(OffMask(i, subscript))));
            rate(logical(~OffMask(i, subscript))) = rate0(subscript(logical(~OffMask(i, subscript))));
            exchange = rand(1, numel(subscript)) < rate;
            OffMask(i, subscript(exchange)) = ~OffMask(i, subscript(exchange));
        end
    end

    %% Real-valued crossover & mutation using fixed mutate rate (1/D)
    if any(Problem.encoding ~= 4)
        OffDec = RealVariation_MutateRateFixed(Problem, ParentDec, OffMask);
        OffDec(:, Problem.encoding == 4) = 1;
    else
        OffDec = ones(N/2, D);
    end
end
