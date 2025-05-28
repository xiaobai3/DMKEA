classdef DMKEA < ALGORITHM
% <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% DMKEA - Multi-stage knowledge-guided evolutionary algorithm
%
% Developed by: Lidan Bai, 2024
%
% Description:
%   DMKEA is a sparse large-scale multi-objective evolutionary algorithm (SLMOEA)
%   designed for problems where only a small subset of decision variables 
%   influence the Pareto optimal front.
%
%   The algorithm combines binary and real-valued variation operators through
%   an adaptive two-stage framework. A dynamic switching probability (sigmoid-based)
%   governs the transition between early-stage variable selection and late-stage value refinement.
%
%   Three knowledge-guided vectors are used to guide the binary search:
%     - Prior vector (pv): estimated via multi-interval sampling during initialization,
%     - Filter guidance vector (fv): captures variance among elite individuals,
%     - Statistical vector (sv): encodes activation frequency of variables across generations.
%
%   Real-valued mutation is restricted to active variables and dynamically adjusts both
%   its mutation rate and distribution index based on the search progress.
%
%   DMKEA is especially effective under extreme sparsity (e.g., θ ≤ 0.01) and was
%   validated on both synthetic and real-world sparse problems. It demonstrates 
%   superior robustness through decoupled treatment of variable selection and 
%   optimization, and avoids hand-crafted stage schedules by employing
%   continuous probability-driven operator selection.
%
% Inspired by:
%   Z. Ding, L. Chen, D. Sun, and X. Zhang, 
%   "A multi-stage knowledge-guided evolutionary algorithm for sparse 
%    multi-objective optimization problems", Swarm and Evolutionary Computation, 
%   2022, 73: 101119.
%
% Note:
%   This algorithm was tested and run within the PlatEMO framework.
%   If using PlatEMO in your research, please cite:
%   Y. Tian, R. Cheng, X. Zhang, and Y. Jin, 
%   "PlatEMO: A MATLAB platform for evolutionary multi-objective optimization",
%   IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87.

    methods
        function main(Algorithm, Problem)
            %% Step 1: Estimate variable importance vector (pv)
            TDec    = [];
            TMask   = [];
            TempPop = [];
            pv      = zeros(1, Problem.D);
            Interval = (Problem.upper - Problem.lower) ./ 4;

            for i = 1 : 1 + 4 * any(Problem.encoding ~= 4)
                for j = 1 : 2
                    Dec = unifrnd(repmat(Problem.lower + Interval*(i-1), Problem.D, 1), ...
                                  repmat(Problem.lower + Interval*i, Problem.D, 1));
                    Dec(:, Problem.encoding == 4) = 1;
                    Mask = eye(Problem.D);
                    Population = Problem.Evaluation(Dec .* Mask);

                    pv = pv + NDSort([Population.objs, Population.cons], inf);

		    % To reduce memory usage and enable parallel execution, we apply incremental
		    % environmental selection after each sampling step. This avoids storing all
		    % temporary populations (TempPop, TDec, TMask) at once, which would otherwise 
		    % require large memory and may crash under high-dimensional settings.
                    [TempPop, TDec, TMask, ~, ~] = SPEA2_EnvironmentalSelection( ...
                        [TempPop, Population], [TDec; Dec], [TMask; Mask], 1000);

                    clear Dec Mask Population;
                end
            end

            pv = pv / (5 * 2);

            %% Step 2: Initialize population with masks guided by pv
            Dec = unifrnd(repmat(Problem.lower, Problem.N, 1), ...
                          repmat(Problem.upper, Problem.N, 1));
            Dec(:, Problem.encoding == 4) = 1;
            Mask = false(Problem.N, Problem.D);
	    
	    % Randomly activate up to 50% of variables for each individual, guided by pv.
	    % This reflects the sparse nature of the problem, where only a small subset of
	    % variables are expected to be relevant. Early-stage sparsity promotes exploration
	    % while avoiding over-activation of redundant variables.
	    % The use of pv ensures that more important variables (lower pv values) are 
	    % more likely to be activated during initialization.
            for i = 1 : Problem.N
                Mask(i, TournamentSelection(2, ceil(0.5 * rand * Problem.D), pv)) = 1;
            end

            Population = Problem.Evaluation(Dec .* Mask);
            [Population, Dec, Mask, FrontNo, CrowdDis] = SPEA2_EnvironmentalSelection( ...
                [Population, TempPop], [Dec; TDec], [Mask; TMask], Problem.N);

            clear TempPop TDec TMask;

            sv = zeros(1, Problem.D);
            Last_temp_num = 0;

            %% Step 3: Main optimization loop
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2, 2 * Problem.N, FrontNo, -CrowdDis);
                First_Mask = Mask(FrontNo == 1, :);
                First_Dec  = Population(FrontNo == 1).decs;

                delta = Problem.FE / Problem.maxFE;

                % Update fv: standard deviation of active decision variables
                fv = std(First_Dec, 0, 1);
                fv(:, Problem.encoding == 4) = sum(First_Mask(:, Problem.encoding == 4), 1);

		% Update sv: statistical vector representing variable activation frequency.
		% A moving average is used, weighted by the number of elite solutions in 
		% the current vs. previous generation, to track sparsity patterns over time.

                temp_num = size(First_Mask, 1);
                temp_vote = sum(First_Mask, 1);
                sv = (Last_temp_num / (Last_temp_num + temp_num)) * sv + ...
                     (temp_num / (Last_temp_num + temp_num)) * (temp_vote / temp_num);
                Last_temp_num = temp_num;
		
		%------------- update pv by sv (optional) ------------------------
		% These lines implement the dynamic update of the prior vector (pv) using 
		% the statistical vector (sv), as described in Equation (6) of the paper. 
		% This update is disabled for synthetic benchmarks (e.g., SMOP1–8) because 
		% pv is reliable from initialization. However, for real-world problems 
		% where variable importance evolves (e.g., Sparse_NN, Sparse_SR), enabling 
		% this update improves guidance quality.
		%
		% Uncomment the line below when applying DMKEA to real-case sparse problems:
		% pv = pv .* (1 - sv) * sqrt(delta) + pv;

                % Determine stage by dynamic probability
                dynamicProb = 1 / (1 + exp(20 * (delta - 0.5)));
	
		% Apply probabilistic switching between two knowledge-guided binary operators:
		% - Operator_PvFv: guided by prior and filter vectors (pv, fv), for early-stage variable identification
		% - Operator_Sv: guided by statistical vector (sv), for later-stage sparsity adaptation
                if rand < dynamicProb
                    [OffDec, OffMask] = Operator_PvFv(Problem, Dec(MatingPool, :), Mask(MatingPool, :), pv, fv, delta);
                else
                    [OffDec, OffMask] = Operator_Sv(Problem, Dec(MatingPool, :), Mask(MatingPool, :), sv);
                end

                Offspring = Problem.Evaluation(OffDec .* OffMask);

                [Population, Dec, Mask, FrontNo, CrowdDis] = SPEA2_EnvironmentalSelection( ...
                    [Population, Offspring], [Dec; OffDec], [Mask; OffMask], Problem.N);
            end
        end
    end
end
