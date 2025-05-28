function [Population, Dec, Mask, FrontNo, CrowdDis] = SPEA2_EnvironmentalSelection(Population, Dec, Mask, N)
% SPEA2_EnvironmentalSelection - Environmental selection used in MSKEA
%
% This function selects N individuals from the current population based on
% non-dominated sorting and truncation, as used in MSKEA.
%
% Inputs:
%   - Population : current population (struct with .objs field)
%   - Dec        : decision variable matrix
%   - Mask       : binary masks for decision variables
%   - N          : number of individuals to select
%
% Outputs:
%   - Population : selected population
%   - Dec        : selected decision variables
%   - Mask       : selected masks
%   - FrontNo    : front number after non-dominated sorting
%   - CrowdDis   : crowding distance
%
% Written by Lei Chen (MSKEA), minor edits by Lidan Bai (2024)

    %% Remove duplicate solutions
    [~, uni] = unique(Population.objs, 'rows');
    Population = Population(uni);
    Dec        = Dec(uni, :);
    Mask       = Mask(uni, :);
    N          = min(N, length(Population));

    %% Non-dominated sorting
    [FrontNo, MaxFNo] = NDSort(Population.objs, N);
    Next = FrontNo < MaxFNo;

    %% Normalize objectives for truncation
    PopObj = Population.objs;
    fmax   = max(PopObj(FrontNo == 1, :), [], 1);
    fmin   = min(PopObj(FrontNo == 1, :), [], 1);
    PopObj = (PopObj - fmin) ./ (fmax - fmin + 1e-12);  % Avoid zero division

    %% Truncation if needed
    Last = find(FrontNo == MaxFNo);
    del  = Truncation(PopObj(Last, :), length(Last) - N + sum(Next));
    Next(Last(~del)) = true;

    %% Final selection
    Population = Population(Next);
    Dec        = Dec(Next, :);
    Mask       = Mask(Next, :);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdingDistance(Population.objs, FrontNo);
end

function Del = Truncation(PopObj, K)
% Truncation - Removes K individuals with minimum crowding

    Distance = pdist2(PopObj, PopObj);
    Distance(logical(eye(size(Distance)))) = inf;
    Del = false(1, size(PopObj, 1));

    while sum(Del) < K
        Remain = find(~Del);
        Temp   = sort(Distance(Remain, Remain), 2);
        [~, idx] = sortrows(Temp);
        Del(Remain(idx(1))) = true;
    end
end
