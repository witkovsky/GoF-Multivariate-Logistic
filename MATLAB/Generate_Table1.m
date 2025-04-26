function [Results, MeansTable, StdsTable] = Generate_Table1(N, dims, sample_sizes, a_vals, kMax, isStandardized)
% GENERATE_TABLE1 Computes empirical expectations (mean ± std) of Rn,a
%   statistics, based on Popović, Mijanović and Witkovský (2025).
%
% SYNTAX:
%   [Results, MeansTable, StdsTable] = Generate_Table1()
%   [Results, MeansTable, StdsTable] = Generate_Table1(N, dims, sample_sizes, a_vals, kMax, isStandardized)
%
% INPUT (optional):
%   N               - number of Monte Carlo repetitions (default: 100)
%   dims            - vector of dimensions (default: [2, 3, 4])
%   sample_sizes    - vector of sample sizes (default: [50, 100])
%   a_vals          - vector of weight parameters (default: [1, 1.5, 2, 2.5, 3, 3.5, 5])
%   kMax            - maximum number of terms in I3 expansion (default: 100)
%   isStandardized  - logical, whether to standardize data (default: true)
%
% OUTPUT:
%   Results         - cell array with means and stds
%   MeansTable      - table with empirical means
%   StdsTable       - table with empirical standard deviations
%
% EXAMPLE:
%   [Results, MeansTable, StdsTable] = Generate_Table1();
%   disp(MeansTable);
%   disp(StdsTable);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025

%% Default settings
if nargin < 1 || isempty(N), N = 100; end
if nargin < 2 || isempty(dims), dims = [2, 3, 4]; end
if nargin < 3 || isempty(sample_sizes), sample_sizes = [50, 100]; end
if nargin < 4 || isempty(a_vals), a_vals = [1, 1.5, 2, 2.5, 3, 3.5, 5]; end
if nargin < 5 || isempty(kMax), kMax = 100; end
if nargin < 6 || isempty(isStandardized), isStandardized = true; end

rng(123)  % Fix random seed for reproducibility

%% Preallocate storage
Results = cell(length(dims), length(sample_sizes), length(a_vals));

%% Display progress
fprintf('Starting simulations...\n')

%% Main Monte Carlo Loop
for id = 1:length(dims)
    d = dims(id);
    for in = 1:length(sample_sizes)
        n = sample_sizes(in);
        for ia = 1:length(a_vals)
            a = a_vals(ia);

            fprintf('Simulating for d = %d, n = %d, a = %.1f...\n', d, n, a)

            Rna_vals = zeros(N, 1);

            for rep = 1:N
                % Generate sample from ML(0, I_d)
                X = randML(n, zeros(1, d), eye(d));

                % Compute Rn,a statistic
                Rna_vals(rep) = Rstat(X, a, kMax, isStandardized);
                % Alternative: Rna_vals(rep) = Rstat_direct(X, a, 5, 10, isStandardized);
            end

            % Store mean ± std
            mean_val = mean(Rna_vals);
            std_val = std(Rna_vals);
            Results{id, in, ia} = [mean_val, std_val];
        end
    end
end

fprintf('Simulation completed!\n')

%% Create Tables for Means and Standard Deviations
headers = {};
mean_values = [];
std_values = [];

for id = 1:length(dims)
    for in = 1:length(sample_sizes)
        colname = sprintf('d%d_n%d', dims(id), sample_sizes(in));
        headers{end+1} = colname; 

        col_mean = zeros(length(a_vals), 1);
        col_std  = zeros(length(a_vals), 1);

        for ia = 1:length(a_vals)
            res = Results{id, in, ia};
            col_mean(ia) = res(1);
            col_std(ia)  = res(2);
        end

        mean_values = [mean_values col_mean]; 
        std_values  = [std_values  col_std];  
    end
end

rowNames = strcat('a', strrep(string(a_vals), '.', '_'));
MeansTable = array2table(mean_values, 'VariableNames', headers, 'RowNames', rowNames);
StdsTable  = array2table(std_values,  'VariableNames', headers, 'RowNames', rowNames);

%% Display
disp('Empirical Means:')
disp(MeansTable)

disp('Empirical Standard Deviations:')
disp(StdsTable)

end