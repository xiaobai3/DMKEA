%% DMKEA Benchmark Runner
% This script runs DMKEA on a set of sparse LSMOP problems using PlatEMO.

%% Add PlatEMO to MATLAB path
addpath('path_to_platemo_folder');  % <-- Replace with actual path

%% Algorithm and problem settings
algorithm = @DMKEA;  % DMKEA main algorithm
problems  = {'SMOP1','SMOP2','SMOP3','SMOP4','SMOP5','SMOP6','SMOP7','SMOP8'};

%% Common parameters
popSize = 100;
D       = 2000;
maxFE   = 300000; % 150D
theta   = 0.01;  % Sparsity

%% Run on each problem
for i = 1:length(problems)
    problem = str2func(problems{i});
    platemo('algorithm', algorithm, ...
            'problem',  problem, ...
            'N',        popSize, ...
            'D',        D, ...
            'maxFE',    maxFE, ...
            'theta',    theta, ...
            'save',     1);  % Set to 1 to save results
end
