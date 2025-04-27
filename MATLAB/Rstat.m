function [Rna, Z] = Rstat(X, a, kMax, isStandardized)
% RSTAT Computes the goodness-of-fit test statistic R_n,a for testing
%  the hypothesis that the X come from a multivariate logistic distribution,
%  based on Theorems 2-4 in Popović, Mijanović and Witkovský (2025).  
%
% The statistic is based on the analytical formula:
%     R_n,a = n * (I1 + I2 - I3),
% where:
%   I1 - empirical integral term
%   I2 - theoretical term under H0
%   I3 - mixed term (empirical vs. theory)
%
% SYNTAX:
%   [R, Z] = Rstat(X, a, kMax, isStandardized)
%
% INPUT:
%   X      - (n x d) matrix of observations (non-standardized)
%   a      - weight function parameter, default is a = 2
%   kMax   - maximum number of terms in series for I3 (default: 100)
%   isStandardized - required standardization of the X (default: true)
%
% OUTPUT:
%   R      - test statistic value
%   Z      - (n x d) matrix of the used (possibly standardized) observations
%
% EXAMPLES:
%   %% Example 1: Bivariate case (standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   kMax = 100;
%   isStandardized = true;
%   [Rna, Z] = Rstat(X, a, kMax, isStandardized);
%   disp(Rna)
%
%   %% Example 2: Bivariate case (non-standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   kMax = 100;
%   isStandardized = false;
%   [Rna, Z] = Rstat(X, a, kMax, isStandardized);
%   disp(Rna)
%
%   %% Example 3: Trivariate case (standardized)
%   n = 100;
%   mu = [0 0 0];
%   Sigma = [1 0.9 0.1; 0.9 1 0.5; 0.1 0.5 1];
%   X = randML(n, mu, Sigma);
%   a = 2;
%   [Rna, Z] = Rstat(X, a);
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
%   [Rna, Z] = Rstat(X, a);
%   disp(Rna)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 27-Apr-2025 09:41:42

%% Input validation
narginchk(1, 4);
if nargin < 2 || isempty(a), a = 2; end
if nargin < 3 || isempty(kMax), kMax = 100; end
if nargin < 4 || isempty(isStandardized), isStandardized = true; end

[n, d] = size(X);

%% Standardize X if isStandardized = true and n > d
if isStandardized
    if n > d
        % Estimate mean and covariance
        mu_hat = mean(X, 1);         % 1 x d vector (row vector of means)
        Sigma_hat = cov(X);          % d x d covariance matrix
        % Cholesky decomposition (upper triangular)
        C = chol(Sigma_hat, 'upper');   % C' * C = Sigma_hat
        % Standardize the sample 
        % Set X such that X' = Sigma_hat^{-1/2}*(X-mu_hat) ]
        % Solving: C * X' = (X - mu_hat)'
        % If standardized rename X
        X = (X - mu_hat) / C;     
    end
end
% Data used (standardized or non-standardized)
Z = X;

%% Compute integrals
I1 = Integral1(Z, a, d, n);
I2 = Integral2(a, d);
I3 = Integral3(Z, a, d, n, kMax);

%% Final test statistic
Rna = n * (I1 + I2 - I3);
end

%% === Integral I1 ===
function I1 = Integral1(X, a, d, n)
% Computes I1 = empirical integral over pairwise distances
if d == 2
    const  = 2 * pi * a / n^2;
    const0 = (2 * pi) / (n * a^2);
    I1 = 0;
    for j = 1:n
        for k = (j+1):n
            r = norm(X(j,:) - X(k,:));
            I1 = I1 + 1 / (r^2 + a^2)^(3/2);
        end
    end
    I1 = const0 + const * 2 * I1;
else
    const = (2 * pi^(d/2) * gamma(d)) / (n^2 * a^d * gamma(d/2));
    I1 = 0;
    for j = 1:n
        for k = (j+1):n
            r2 = norm(X(j,:) - X(k,:))^2;
            I1 = I1 + Hypergeom2F1(d/2, d/2 + 0.5, d/2, -r2 / a^2);
        end
    end
    I1 = const * (n + 2*I1);
end
end

%% === Integral I2 ===
function [I2, RINF] = Integral2(a, d)
% Computes I2 = theoretical expectation under H0 and RINF
z = a / (2 * sqrt(3));
const = (2^(1 - d) * pi^(d/2) * gamma(d + 2)) / (3^(d/2) * gamma(d/2));
const2 = (2 * pi^(d/2) * gamma(d)) / (a^d * gamma(d/2));

I2 = hurwitzZeta(d + 1, z) - z * hurwitzZeta(d + 2, z);
I2 = const * I2;

RINF = const2 - I2;
end

%% === Integral I3 ===
function I3 = Integral3(X, a, d, n, kMax)
% Computes I3 = integral involving X vs. CF (series over k)
normData2 = sum(X.^2, 2);
sqrt3 = sqrt(3);
I3 = zeros(n, 1);

for k = 0:kMax
    denom = (a + sqrt3 * (2 * k + 1))^2;
    hyp = Hypergeom2F1(d/2 + 0.5, d/2 + 1, d/2, -normData2 / denom);
    I3 = I3 + hyp ./ denom.^(d/2 + 0.5);
end

const = (8 * sqrt3 * pi^(d/2) * gamma(d + 1)) / (n * gamma(d/2));
I3 = const * sum(I3);
end