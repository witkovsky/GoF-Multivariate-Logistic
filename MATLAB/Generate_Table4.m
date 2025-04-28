function [EmpiricalQuantiles, EmpiricalTables] = ...
    Generate_Table4(N, dims, a_vals, n_samples, alpha, kMax, Tmax, grid_pts, isStandardized)
% GENERATE_TABLE4 Computes empirical quantiles of R_{n,a} statistics for given
%   significance levels. Fast computation: uses Rstat for small n, Rstat_direct
%   for large n. 
%
% SYNTAX:
%   [EmpiricalQuantiles, EmpiricalTables] = Generate_Table4()
%
%   [EmpiricalQuantiles, EmpiricalTables] = Generate_Table4(N, dims, ...
%          a_vals, n_samples, alpha, kMax, Tmax, grid_pts, isStandardized)
% %
% INPUT (optional):
%   N               - number of MC repetitions for all n (default: 1000)
%   dims            - vector of dimensions (default: [2, 3, 4])
%   a_vals          - vector of weight parameters (default: [1, 1.5, 2, 2.5, 3, 3.5, 5])
%   n_samples       - vector of sample sizes (default: [30, 50, 75, 100, 150, 200, 300, 10000])
%   alpha           - vector of significance levels (default: 0.05)
%   kMax            - maximum terms for expansion in Rstat (default: 100)
%   Tmax            - maximum integration limit for Rstat_direct (default: 5)
%   grid_pts        - number of grid points for Rstat_direct (default: 50)
%   isStandardized  - whether to standardize data before test (default: true)
%
% OUTPUT:
%   EmpiricalQuantiles - structure with 3D arrays of empirical quantiles (a x n x alpha)
%   EmpiricalTables    - structure with tables for each dimension and alpha level
%
%
% EXAMPLE:
%   [EmpiricalQuantiles, EmpiricalTables] = Generate_Table4();
%   disp(EmpiricalQuantiles);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 27-Apr-2025

%% Default settings
if nargin < 1 || isempty(N), N = 1000; end
if nargin < 2 || isempty(dims), dims = [2, 3, 4]; end
if nargin < 3 || isempty(a_vals), a_vals = [1, 1.5, 2, 2.5, 3, 3.5, 5]; end
if nargin < 4 || isempty(n_samples), n_samples = [30, 50, 75, 100, 150, 200, 300, 10000]; end
if nargin < 5 || isempty(alpha), alpha = 0.05; end
if nargin < 6 || isempty(kMax), kMax = 100; end
if nargin < 7 || isempty(Tmax), Tmax = 5; end
if nargin < 8 || isempty(grid_pts), grid_pts = 50; end
if nargin < 9 || isempty(isStandardized), isStandardized = true; end

rng(123); % Fix random seed for reproducibility

%% Preallocate
EmpiricalQuantiles = struct();
EmpiricalTables = struct();

%% Main Simulation Loop
fprintf('Starting simulations for empirical quantiles...\n')

for id = 1:length(dims)
    d = dims(id);
    QuantileMatrix = zeros(length(a_vals), length(n_samples), length(alpha));

    for ia = 1:length(a_vals)
        a = a_vals(ia);

        for in = 1:length(n_samples)
            n = n_samples(in);
            fprintf('Dimension d = %d, weight a = %.1f, sample size n = %d\n', d, a, n)

            % Adjust grid points for large n
            if n >= 10000
                current_grid_pts = max(5, round(grid_pts/d));
            elseif n >= 1000
                current_grid_pts = max(10, round(grid_pts/d));
            else
                current_grid_pts = grid_pts;
            end

            Rna_vals = zeros(N, 1);

            for rep = 1:N
                X = randML(n, zeros(1,d), eye(d));
                if n <= 100
                    Rna_vals(rep) = Rstat(X, a, kMax, isStandardized);
                else
                    Rna_vals(rep) = Rstat_direct(X, a, Tmax, ...
                        current_grid_pts, isStandardized);
                end
            end

            % Compute empirical quantiles for all alpha levels
            for ialpha = 1:length(alpha)
                QuantileMatrix(ia, in, ialpha) = ...
                    quantile(Rna_vals, 1 - alpha(ialpha));
            end
        end
    end

    EmpiricalQuantiles.(sprintf('d%d',d)) = QuantileMatrix;
    save Table4

    % Create tables for each alpha level
    for ialpha = 1:length(alpha)
        rowNames = strcat('a', strrep(string(a_vals), '.', '_'));
        colNames = strcat('n', string(n_samples));
        MatrixData = QuantileMatrix(:,:,ialpha);
        EmpiricalTables.(sprintf('d%d_alpha%d',d,1000*alpha(ialpha))) = ...
            array2table(MatrixData, 'VariableNames', colNames, 'RowNames', rowNames);
    end

    save Table4
end

%% Display
for id = 1:length(dims)
    d = dims(id);
    for ialpha = 1:length(alpha)
        fprintf('\nResults for dimension d = %d, significance level alpha = %.3f:\n', d, alpha(ialpha));
        disp(EmpiricalTables.(sprintf('d%d_alpha%d',d,1000*alpha(ialpha))));
    end
end

end
