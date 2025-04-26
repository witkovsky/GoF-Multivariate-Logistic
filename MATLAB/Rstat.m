function [R, RINF] = Rstat(data, a, kMax, isStandardized)
% RSTAT Computes the goodness-of-fit test statistic R_n,a for testing
% the hypothesis that the data come from a multivariate logistic distribution.
%
% The statistic is based on the analytical formula:
%     R_n,a = n * (I1 + I2 - I3),
% where:
%   I1 - empirical integral term
%   I2 - theoretical term under H0
%   I3 - mixed term (data vs. theory)
%
% SYNTAX:
%   [R, RINF] = Rstat(data, a, kMax, isStandardized)
%
% INPUT:
%   data   - (n x d) matrix of observations (standardized)
%   a      - weight function parameter, default is a = 2
%   kMax   - maximum number of terms in series for I3 (default: 100)
%   isStandardized - standardization of the data (default: true)
%
% OUTPUT:
%   R      - test statistic value
%   RINF   - expected value under H0 (used for calibration)
%
% EXAMPLES:
%   %% Example 1: Bivariate case (standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   data = randML(n, mu, Sigma);
%   a = 2;
%   kMax = 100;
%   isStandardized = true;
%   R = Rstat(data, a, kMax, isStandardized)
%
%   %% Example 2: Bivariate case (non-standardized)
%   rng(1)
%   n = 100;
%   mu = [0 0];
%   Sigma = [1 0.9; 0.9 2];
%   data = randML(n, mu, Sigma);
%   a = 2;
%   kMax = 100;
%   isStandardized = false;
%   R = Rstat(data, a, kMax, isStandardized)
%
%   %% Example 3: Trivariate case (standardized)
%   n = 100;
%   mu = [0 0 0];
%   Sigma = [1 0.9 0.1; 0.9 1 0.5; 0.1 0.5 1];
%   data = randML(n, mu, Sigma);
%   a = 2;
%   R = Rstat(data, a)
%
%   %% Example 4: 4-dimensional case (standardized)
%   n = 100;
%   mu = [1 2 3 4];
%   Sigma = [1 0.9 0.1 0;
%            0.9 1 0.5 0;
%            0.1 0.5 1 0;
%            0 0 0 2];
%   data = randML(n, mu, Sigma);
%   a = 2;
%   R = Rstat(data, a)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2025 09:41:42

%% Input validation
narginchk(1, 4);
if nargin < 2 || isempty(a), a = 2; end
if nargin < 3 || isempty(kMax), kMax = 100; end
if nargin < 4 || isempty(isStandardized), isStandardized = true; end

[n, d] = size(data);

%% Standardize data if isStandardized = true
if isStandardized
    muHat = mean(data);
    SigmaHat = cov(data);
    data = (sqrtm(inv(SigmaHat)) * (data - muHat)')';
end

%% Compute integrals
I1 = Integral1(data, a, d, n);
[I2, RINF] = Integral2(a, d);
I3 = Integral3(data, a, d, n, kMax);

%% Final test statistic
R = n * (I1 + I2 - I3);
end

%% === Integral I1 ===
function I1 = Integral1(data, a, d, n)
% Computes I1 = empirical integral over pairwise distances
if d == 2
    const = 2 * pi * a / n^2;
    I1 = 0;
    for j = 1:n
        for k = 1:n
            r = norm(data(j,:) - data(k,:));
            I1 = I1 + 1 / (r^2 + a^2)^(3/2);
        end
    end
    I1 = const * I1;
else
    const = (2 * pi^(d/2) * gamma(d)) / (n^2 * a^d * gamma(d/2));
    I1 = 0;
    for j = 1:n
        for k = 1:n
            r2 = norm(data(j,:) - data(k,:))^2;
            I1 = I1 + Hypergeom2F1(d/2, d/2 + 0.5, d/2, -r2 / a^2);
        end
    end
    I1 = const * I1;
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
function I3 = Integral3(data, a, d, n, kMax)
% Computes I3 = integral involving data vs. CF (series over k)
normData2 = sum(data.^2, 2);
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