% TABLE 1
% SCRIPT: Generate values for Table 1 - Empirical expectations of Rn,a.
% Based on Popović, Mijanović and Witkovský (2025).
%

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025 09:41:42

clear; 
close all; 
clc

%% Settings
N = 100;                               % Number of Monte Carlo repetitions
dims = [2, 3, 4];                      % Dimensions
sample_sizes = [50, 100];              % Sample sizes
a_vals = [1, 1.5, 2, 2.5, 3, 3.5, 5];  % Values of a
kMax = 100;                            % maximum number of terms in I3 
isStandardized = true;                 % standardization of the data
rng(123)                               % Reset random number generator 
                                       % for reproducibility


% Preallocate results
Results = cell(length(dims), length(sample_sizes), length(a_vals));

% Display progress
fprintf('Starting simulations...\n')

%% Main Loop
for id = 1:length(dims)
    d = dims(id);
    for in = 1:length(sample_sizes)
        n = sample_sizes(in);
        for ia = 1:length(a_vals)
            a = a_vals(ia);
            
            Rna_vals = zeros(N,1); % store Rna values
            
            fprintf('Simulating for d = %d, n = %d, a = %.1f...\n', d, n, a)
            
            for rep = 1:N
                % Generate sample from ML(0,I)
                X = randML(n, zeros(1,d), eye(d)); % you must have randML implemented
                
                % Compute Rn,a statistic
                Rna_vals(rep) = Rstat(X, a, kMax, isStandardized);
                % Alternatively, use:
                % Rna_vals(rep) = Rstat_direct(X, a, 5, 10, standardize);
            end
            
            % Store mean ± std
            mean_val = mean(Rna_vals);
            std_val = std(Rna_vals);
            Results{id, in, ia} = [mean_val, std_val];
        end
    end
end

fprintf('Simulation completed!\n')

%% Display Results in a Table Format
disp('Empirical Means ± Standard Deviations:')
for id = 1:length(dims)
    for in = 1:length(sample_sizes)
        fprintf('Dimension d = %d, Sample size n = %d\n', dims(id), sample_sizes(in));
        for ia = 1:length(a_vals)
            res = Results{id, in, ia};
            fprintf('a = %.1f: Mean = %.5f, Std = %.5f\n', a_vals(ia), res(1), res(2));
        end
        fprintf('\n');
    end
end
