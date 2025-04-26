function [Rna, Z] = Rstat_direct(X, a, Tmax, grid_pts, isStandardized)
% RSTAT_DIRECT Computes the goodness-of-fit test statistic R_n,a
%   for testing the hypothesis that X comes from a multivariate logistic
%   distribution. The statistic is based on direct numerical integration
%   over a grid, see Eq. (8) in Popović, Mijanović and Witkovský (2025). 
%
% SYNTAX:
%   [R, Z] = Rstat_direct(X, a, Tmax, grid_pts, isStandardized)
%
% INPUT:
%   X              - (n x d) matrix of observations (non-standardized)
%   a              - weight function parameter, default is a = 3
%   Tmax           - range parameter for grid [-Tmax, Tmax] (default: 5)
%   grid_pts       - number of grid points per dimension (default: 100)
%   isStandardized - logical, whether to standardize X (default: true)
%
% OUTPUT:
%   R              - test statistic value
%   Z              - (n x d) matrix of the used (possibly standardized)
%                    observations 
%
% DEPENDENCIES:
%   emp_cf.m, cf_logistic.m
%
% EXAMPLES:
%   %% Example 1: Bivariate case (standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   Tmax = 5;
%   grid_pts = 100;
%   isStandardized = true;
%   [Rna, Z] = Rstat_direct(X, a, Tmax, grid_pts, isStandardized);
%   disp(Rna)
%
%   %% Example 1: Bivariate case (non-standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   Tmax = 5;
%   grid_pts = 100;
%   isStandardized = false;
%   [Rna, Z] = Rstat_direct(X, a, Tmax, grid_pts, isStandardized);
%   disp(Rna)
%
%   %% Example 3: Trivariate case (non-standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0 0];
%   Sigma = [1 0.9 0.1; 0.9 1 0.5; 0.1 0.5 1];
%   X = randML(n, mu, Sigma);
%   [Rna, Z] = Rstat_direct(X)
%   disp(Rna)
%
%   %% Example 4: 4-dimensional case (standardized)
%   n = 100;
%   mu = [1 2 3 4];
%   Sigma = [1 0.9 0.1 0;
%            0.9 1 0.5 0;
%            0.1 0.5 1 0;
%            0 0 0 2];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   Tmax = 5;
%   grid_pts = 10;
%   isStandardized = true;
%   [Rna, Z] = Rstat_direct(X, a, Tmax, grid_pts, isStandardized);
%   disp(Rna)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025

%% Input validation
narginchk(1, 5);
if nargin < 2 || isempty(a), a = 2; end
if nargin < 3 || isempty(Tmax), Tmax = 5; end
if nargin < 4 || isempty(grid_pts), grid_pts = 100; end
if nargin < 5 || isempty(isStandardized), isStandardized = true; end

[n, d] = size(X);

%% Standardize X if isStandardized = true and n > d
if isStandardized
    if n > d
        mu_hat = mean(X, 1);
        Sigma_hat = cov(X);
        C = chol(Sigma_hat, 'upper');
        X = (X - mu_hat) / C;
    end
end
Z = X;

%% Create grid of t-points
vec = linspace(-Tmax, Tmax, grid_pts);
grid_cell = cell(1, d);
[grid_cell{:}] = ndgrid(vec);
t_mat = reshape(cat(d+1, grid_cell{:}), [], d);

%% Compute empirical and theoretical characteristic functions
phi_n = emp_cf(Z, t_mat);
phi_0 = cf_logistic(t_mat);

%% Compute squared differences
diff_sq = (real(phi_n) - phi_0).^2 + (imag(phi_n)).^2;

%% Compute weights
weights = exp(-a * sqrt(sum(t_mat.^2, 2)));

%% Volume element for integration
volume_element = (2 * Tmax / grid_pts)^d;

%% Final R_n,a computation
Rna = n * sum(diff_sq .* weights) * volume_element;

end
%% function phi_n = emp_cf(Z, t)
function phi_n = emp_cf(Z, t)
% EMP_CF Computes the empirical characteristic function (ECF) of sample Z evaluated at points t.
%
% INPUT:
%   Z - n x d matrix, sample data
%   t - m x d matrix, evaluation points (each row is a t-vector)
%
% OUTPUT:
%   phi_n - m x 1 complex vector, empirical CF evaluated at each t
%
% EXAMPLE:
%   Z = randn(100,2);
%   t = randn(50,2);
%   phi_n = emp_cf(Z, t);

[~, d] = size(Z);
[~, d_check] = size(t);

if d ~= d_check
    error('Dimension mismatch between sample Z and evaluation points t.');
end

% Compute t * Z' inner products
% Resulting in an m x n matrix: each element (i,j) is t_i' * Z_j
inner_products = t * Z';

% Exponentiate and average over the sample
phi_n = mean(exp(1i * inner_products), 2);  % mean over columns (i.e., over n samples)

end
%% function res = cf_logistic(t)
function res = cf_logistic(t)
% CF_LOGISTIC Computes the characteristic function of the multivariate logistic distribution.
%
% INPUT:
%   t - n x d matrix, where each row is a vector at which the CF is evaluated
%
% OUTPUT:
%   res - n x 1 vector, values of the characteristic function at each t
%
% EXAMPLE:
%   t = randn(100, 2);
%   res = cf_logistic(t);

%% Compute norms
norm_t = sqrt(sum(t.^2, 2));  % Euclidean norm for each row

%% Compute characteristic function
sqrt3 = sqrt(3);
res = sqrt3 * norm_t ./ sinh(sqrt3 * norm_t);

%% Handle 0/0 at norm_t = 0 (limit is 1)
res(isnan(res)) = 1;

end

