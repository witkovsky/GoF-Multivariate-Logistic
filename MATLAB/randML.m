function [sample, Y, U, Ucorr, correction_coef, mu, Sigma] = randML(n, mu, Sigma)
% randML generates a random sample of size n from the multivariate logistic
% distribution ML_d(mu, Sigma), as defined by Popović, Mijanović, and Witkovský (2024).
%
% The generation follows the stochastic representation:
%     X = (E[U^2])^{-1/2} * U * Y + mu,
% where:
%     - Y ~ N(0, Sigma) (d-dimensional multivariate normal),
%     - U ~ Kolmogorov–Smirnov-type distribution (Eq. 1 in the paper),
%     - mu is a location vector in R^d,
%     - Sigma is a symmetric positive definite d x d scatter matrix.
%
% The resulting distribution has the characteristic function:
%     φ_X(t) = exp(i * mu' * t) * sqrt(3 * t' * Sigma * t) / sinh(sqrt(3 * t' * Sigma * t))
%
% SYNTAX:
%   [sample, Y, U, Ucorr, correction_coef, mu, Sigma] = randML(n, mu, Sigma)
%
% INPUTS:
%   n      - Number of random samples
%   mu     - Location vector (1 × d)
%   Sigma  - Scatter matrix (d × d, positive definite)
%
% OUTPUTS:
%   sample           - n × d matrix of multivariate logistic random vectors
%   Y                - n × d matrix of standard normal vectors
%   U                - n × 1 vector of Kolmogorov-Smirnov-type variables
%   Ucorr            - Scaled U to ensure unit variance (for standardization)
%   correction_coef  - Scaling constant (E[U^2])^{-1/2}
%   mu               - Re-used location vector (after defaulting if needed)
%   Sigma            - Re-used scatter matrix (after defaulting if needed)
%
% EXAMPLES:
%   % Example: d = 2
%   n = 100; mu = [1, 2]; Sigma = [1 0.75; 0.75 5];
%   X = randML(n, mu, Sigma);
%   mean(X), cov(X)
%
%   % Example: d = 3
%   mu = [1, 2, 3];
%   Sigma = [1 0.75 -0.25; 0.75 2 0.1; -0.25 0.1 3];
%   X = randML(100, mu, Sigma);
%   mean(X), cov(X)

% (c) Viktor Witkovsky, edited and documented for reproducibility
% Ver.: 'Updated on 2025-04-16'

%% Input checking and defaults
narginchk(1, 3);
if nargin < 2 || isempty(mu), mu = [0, 0]; end
if nargin < 3 || isempty(Sigma), Sigma = eye(length(mu)); end

%% Step 1: Generate standard multivariate normal vectors
% Y ~ N(0, Sigma), d-dimensional
Y = mvnrnd(zeros(1, length(mu)), Sigma, n);

%% Step 2: Generate Kolmogorov–Smirnov random variable U
% U ~ Kolmogorov–Smirnov-type distribution as in Eq. (1) of the paper
U = KSinv(rand(n, 1));

%% Step 3: Normalize U to unit variance using correction coefficient
% E[U^2] = pi^2 / (3 * 4), correction = 1 / sqrt(E[U^2])
correction_coef = 1 / sqrt(pi^2 / (3 * 4));
Ucorr = correction_coef * U;

%% Step 4: Construct sample
% Each row of sample: mu + Ucorr(i) * Y(i,:)
sample = bsxfun(@plus, Ucorr .* Y, mu);

end

%% Function KSinv
function x = KSinv(p)
% KSinv Numerically computes the inverse CDF of the Kolmogorov-type distribution
% Defined via the series in Eq. (1) of the paper using root-finding (fzero)

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
%% Function KScdf
function p = KScdf(x)
% KScdf Computes the Kolmogorov-type cumulative distribution function
% Approximation via the infinite series:
%   P(U ≤ u) = sqrt(2*pi) * sum_{k=1}^∞ exp(-((2k−1)^2 * π^2) / (8u^2)) / u

if nargin < 1
    error('KScdf requires at least one input argument.');
end

p = zeros(size(x));
idx = find(x > 0);
xpos = x(idx);

pval = zeros(size(xpos));
for k = 1:1000
    term = exp(-(2 * k - 1)^2 * pi^2 ./ (8 * xpos.^2)) ./ xpos;
    pval = pval + sqrt(2 * pi) * term;
    if all(abs(term) < 1e-12), break; end
end

p(idx) = pval;
end