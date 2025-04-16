function [sample, Y, U, Ucorr, correction_coef, mu, Sigma] = randML(n, mu, Sigma)
% randML generates a random sample of size n from a multivariate Logistic
% distribution with location vector `mu` and scatter matrix `Sigma`.
%
% SYNTAX:
%   [sample, Y, U, Ucorr, correction_coef] = randML(n, mu, Sigma)
%
% INPUTS:
%   n      - Number of random samples
%   mu     - Mean (location) vector of the distribution
%   Sigma  - Scatter (covariance-like) matrix
%
% OUTPUTS:
%   sample           - n-by-d matrix of multivariate logistic random samples
%   Y                - Standardized normal random vectors
%   U                - Scalar logistic-like random variables (KS-transformed)
%   Ucorr            - Rescaled version of U for variance adjustment
%   correction_coef  - Coefficient used to adjust U
%   mu               - Final mean vector used
%   Sigma            - Final scatter matrix used
%
% EXAMPLES:
%   % Bivariate case
%   n = 100;
%   mu = [1 2];
%   Sigma = [1 0.75; 0.75 5];
%   data = randML(n, mu, Sigma);
%   mean(data), cov(data)
%
%   % Trivariate case
%   n = 100;
%   mu = [1 2 3];
%   Sigma = [1 0.75 -0.25; 0.75 2 0.1; -0.25 0.1 3];
%   data = randML(n, mu, Sigma);
%   mean(data), cov(data)
%
% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '21-Apr-2024 12:02:25'

%% Input validation and defaults
narginchk(1, 3);
if nargin < 2 || isempty(mu), mu = [0, 0]; end
if nargin < 3 || isempty(Sigma), Sigma = eye(length(mu)); end

%% Generate components
% Y ~ N(0, Sigma), d-dimensional multivariate normal
Y = mvnrnd(zeros(1, length(mu)), Sigma, n);

% U ~ Kolmogorov-type distribution, scalar per sample
U = KSinv(rand(n, 1));

% Correction coefficient from Balakrishnan's formulation:
% ensures variance of U matches logistic scale factor
correction_coef = 1 / sqrt(pi^2 / (3 * 4));
Ucorr = correction_coef * U;

% Combine the components: scale each Y(i,:) by Ucorr(i), then shift by mu
sample = bsxfun(@plus, Ucorr .* Y, mu);

end

%% Inverse of Kolmogorov CDF (numerically via root finding)
function x = KSinv(p)
% KSinv Inverse of the Kolmogorov CDF using numerical root-finding

x = zeros(size(p));
x(p < 0 | p > 1) = NaN;
x(p == 0) = 0;
x(p == 1) = Inf;

fun = @(x, pval) KScdf(x) - pval;
options = optimset('TolX', sqrt(eps));

for i = 1:numel(p)
    if p(i) > 0 && p(i) < 1
        x(i) = fzero(@(x) fun(x, p(i)), 0, options);
    end
end
end

%% Kolmogorov CDF approximation (based on infinite series)
function p = KScdf(x)
% KScdf Kolmogorov cumulative distribution function.
% This is a series expansion approximation based on the definition
% commonly found in mathematical literature (e.g., Wikipedia).

if nargin < 1
    error('KScdf requires at least one input argument.');
end

p = zeros(size(x));
idx = find(x > 0);
xpos = x(idx);

% Approximation by series expansion
pval = zeros(size(xpos));
for k = 1:1000
    % Wikipedia-based series
    add = exp(-(2 * k - 1)^2 * pi^2 ./ (8 * xpos.^2)) ./ xpos;
    pval = pval + sqrt(2 * pi) * add;
    if all(abs(add) < 1e-12), break; end
end

p(idx) = pval;
end
