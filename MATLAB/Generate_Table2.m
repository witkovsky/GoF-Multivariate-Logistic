function [results, resultsTable] = Generate_Table2(dims, a_vals, method)
% GENERATE_TABLE2 Computes theoretical expectations E[R_inf,a] for Table 2,
% based on Popović, Mijanović and Witkovský (2025).
% 
% Theoretical values are computed using either:
%   - Rstat_ExpectedValue_Expansion (preferred, fast)
%   - Rstat_ExpectedValue_Integral (alternative, numerical integration)
%
% SYNTAX:
%   results = Generate_Table2()
%   [results, resultsTable] = Generate_Table2(dims, a_vals, method)
%
% INPUT:
%   dims    - vector of dimensions (default: [2, 3, 4])
%   a_vals  - vector of weight parameters (default: [1, 1.5, 2, 2.5, 3, 3.5, 5])
%   method  - 'expansion' (default) or 'integral'
%
% OUTPUT:
%   results      - matrix of computed E[R_inf,a] values
%   resultsTable - table with a values in rows and dimensions in columns
%
% EXAMPLE:
%   [results, tbl] = Generate_Table2();
%   disp(tbl);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025

%% Default settings
if nargin < 1 || isempty(dims)
    dims = [2, 3, 4];
end
if nargin < 2 || isempty(a_vals)
    a_vals = [1, 1.5, 2, 2.5, 3, 3.5, 5];
end
if nargin < 3 || isempty(method)
    method = 'expansion';
end

% Preallocate results
results = zeros(length(a_vals), length(dims));

%% Main computation
for id = 1:length(dims)
    d = dims(id);
    for ia = 1:length(a_vals)
        a = a_vals(ia);

        switch lower(method)
            case 'expansion'
                % Compute using expansion method
                Rinf = Rstat_ExpectedValue_Expansion(d, a);
            case 'integral'
                % Compute using integral method
                Rinf = Rstat_ExpectedValue_Integral(d, a);
            otherwise
                error('Unknown method. Choose "expansion" or "integral".');
        end

        results(ia, id) = Rinf;
    end
end

%% Create table output
resultsTable = array2table(results, 'VariableNames', strcat('d', string(dims)), 'RowNames', strcat('a', strrep(string(a_vals), '.', '_')));

%% Display results
fprintf('Theoretical Expectations E[R_inf,a]:\n');
for id = 1:length(dims)
    fprintf('Dimension d = %d\n', dims(id));
    for ia = 1:length(a_vals)
        fprintf('a = %.1f: E[R_inf,a] = %.6f\n', a_vals(ia), results(ia, id));
    end
    fprintf('\n');
end
end