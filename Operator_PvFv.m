function [OffDec, OffMask] = Operator_PvFv(Problem, ParentDec, ParentMask, pv, fv, delta)
% Operator_PvFv - Binary variation guided by prior (pv) and filter (fv) vectors
%
% This operator is part of the DMKEA framework for sparse multi-objective optimization.
% It uses a combination of pv-guided crossover and a probabilistic mutation strategy,
% dynamically switching between pv- and fv-based guidance depending on the search progress (delta).
%
% Developed by: Lidan Bai, 2024

    %% Parameter setting
    [N, D] = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2, :);
    Parent2Mask = ParentMask(N/2+1:end, :);

    %% Binary crossover guided by pv
    % For each pair, flip one differing bit based on pv: remove high-pv or add low-pv variable
    OffMask = Parent1Mask;
    for i = 1:N/2
        if rand < 0.5
            % Attempt to remove less useful variables
            index = find(Parent1Mask(i, :) & ~Parent2Mask(i, :));
            index = index(TS(-pv(index)));  % Prefer high-pv
            OffMask(i, index) = 0;
        else
            % Attempt to add promising variables
            index = find(~Parent1Mask(i, :) & Parent2Mask(i, :));
            index = index(TS(pv(index)));   % Prefer low-pv
            OffMask(i, index) = Parent2Mask(i, index);
        end
    end

    %% Binary mutation guided by pv or fv
    % Use delta to control the guidance source:
    % - if rand < delta → use fv (filter-guided)
    % - else           → use pv (prior-guided)
    f_vector = double(fv > 0);  % Auxiliary binary vector from fv
    for i = 1:N/2
        if rand < delta
            index = find(OffMask(i, :) ~= f_vector);
            if rand < 0.5
                index = index(TS(-fv(index)));  % Flip bit with high fv → 1
                OffMask(i, index) = 1;
            else
                index = index(TS(fv(index)));   % Flip bit with low fv → 0
                OffMask(i, index) = 0;
            end
        else
            if rand < 0.5
                index = find(OffMask(i, :) & pv > median(pv));
                index = index(TS(-pv(index)));  % Flip bit with high pv → 0
                OffMask(i, index) = 0;
            else
                index = find(~OffMask(i, :) & pv < median(pv));
                index = index(TS(pv(index)));   % Flip bit with low pv → 1
                OffMask(i, index) = 1;
            end
        end
    end

    %% Real-valued variation (if encoding is not pure binary)
    if any(Problem.encoding ~= 4)
        OffDec = RealVariation_MutateRateByActive(Problem, ParentDec, OffMask);  % Your custom operator
        OffDec(:, Problem.encoding == 4) = 1;
    else
        OffDec = ones(N/2, D);
    end
end

function index = TS(pv)
% TS - Binary tournament selection for 1 winner
% Selects index with better value (lower for minimization)

    if isempty(pv)
        index = [];
    else
        index = TournamentSelection(2, 1, pv);
    end
end
