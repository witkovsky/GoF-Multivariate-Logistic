function [CriticalValues, CriticalTable] = Generate_Table3(N, dims, sample_sizes, a_vals, kMax, isStandardized)
% GENERATE_TABLE3 Computes 5% critical values of Rn,a statistics under ML(0,I).
% Based on Popović, Mijanović and Witkovský (2025).
%
% SYNTAX:
%   [CriticalValues, CriticalTable] = Generate_Table3()
%   [CriticalValues, CriticalTable] = Generate_Table3(N, dims, sample_sizes, a_vals, kMax, isStandardized)
%
% INPUT (optional):
%   N               - number of Monte Carlo repetitions (default: 1000)
%   dims            - vector of dimensions (default: [2, 3, 4])
%   sample_sizes    - vector of sample sizes (default: [50, 100])
%   a_vals          - vector of weight parameters (default: [1, 1.5, 2, 2.5, 3, 3.5, 5])
%   kMax            - maximum number of terms in I3 expansion (default: 100)
%   isStandardized  - logical, whether to standardize data (default: true)
%
% OUTPUT:
%   CriticalValues  - cell array with critical values
%   CriticalTable   - table of critical values (5% level)
%
% EXAMPLE:
%   [CriticalValues, CriticalTable] = Generate_Table3();
%   disp(CriticalTable);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025

%% Default settings
if nargin < 1 || isempty(N), N = 1000; end
if nargin < 2 || isempty(dims), dims = [2, 3, 4]; end
if nargin < 3 || isempty(sample_sizes), sample_sizes = [50, 100]; end
if nargin < 4 || isempty(a_vals), a_vals = [1, 1.5, 2, 2.5, 3, 3.5 5]; end
if nargin < 5 || isempty(kMax), kMax = 100; end
if nargin < 6 || isempty(isStandardized), isStandardized = true; end

rng(123)  % Fix random seed for reproducibility

%% Preallocate storage
CriticalValues = cell(length(dims), length(sample_sizes), length(a_vals));

%% Display progress
fprintf('Starting critical value simulations...\n')

%% Main Loop
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

            % Compute empirical 95% quantile (critical value)
            critical_val = quantile(Rna_vals, 0.95);

            CriticalValues{id, in, ia} = critical_val;
        end
    end
end

fprintf('Simulation completed!\n')

%% Create Table of Critical Values
headers = {};
critical_values = [];

for id = 1:length(dims)
    for in = 1:length(sample_sizes)
        colname = sprintf('d%d_n%d', dims(id), sample_sizes(in));
        headers{end+1} = colname; %#ok<AGROW>

        col_crit = zeros(length(a_vals), 1);
        for ia = 1:length(a_vals)
            val = CriticalValues{id, in, ia};
            col_crit(ia) = val;
        end

        critical_values = [critical_values col_crit]; %#ok<AGROW>
    end
end

rowNames = strcat('a', strrep(string(a_vals), '.', '_'));
CriticalTable = array2table(critical_values, 'VariableNames', headers, 'RowNames', rowNames);

%% Display
disp('Empirical 5% Critical Values:')
disp(CriticalTable)

end