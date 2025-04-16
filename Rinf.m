function [R, R0, terms] = Rinf(a, d, N, M)
% RINF Computes the expected value E[R_inf,a] as given in Theorem 8
% for the goodness-of-fit test statistic based on the multivariate logistic
% distribution. 
%
% SYNTAX:
%   [R, R0, terms] = Rinf(a, d, N, M)
%
% INPUT:
%   a - positive scalar (weight parameter)
%   d - positive integer (dimension)
%   N - number of explicit terms in sums S1, S2, S3 (default: 1000)
%   M - number of Bernoulli correction terms in S1, S2, S3 (default: 2)
%
% OUTPUT:
%   R - value of E[R_inf,a] based on Theorem 8 expression
%
% REQUIRES:
%   HurwitzZeta.m, SumS1.m, SumS2.m, SumS3.m
%
% EXAMPLE:
%   Rinf(2, 3)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '16-Apr-2025 16:03:44'

%% ALGORITHM
%% ALGORITHM
if nargin < 3, N = 1000; end
if nargin < 4, M = 2; end

% Validate inputs
if a <= 0 || d <= 0 || mod(d,1) ~= 0
    error('Input a must be > 0 and d must be a positive integer.');
end

% Constants
sqrt3 = sqrt(3);
z = a / (2 * sqrt3);
pi_d2 = pi^(d / 2);
gamma_d2 = gamma(d / 2);

% Hurwitz zeta terms
hz1 = HurwitzZeta(d + 1, z);
hz2 = HurwitzZeta(d + 2, z);

% First two terms (standard part)
term1 = (2 * pi_d2 * gamma(d)) / (gamma_d2 * a^d);
term2 = - (2^(1 - d) * pi_d2 * gamma(d + 2)) / (3^(d / 2) * gamma_d2) * (hz1 - z * hz2);

% Prepare indices for S1, S2, S3
% m_vals = [d-1, d, d+1, d+2, d+3];

% Evaluate sums with selected parameters N and M
S1 = @(m) SumS1(m, a, N, M);
S2 = @(m) SumS2(m, a, N, M);
S3 = @(m) SumS3(m, a, N, M);

% Additional terms in the full expansion
term3  = (48 * pi_d2 * gamma(d+2)) / (gamma_d2 * 3^((d+2)/2)) * (S1(d+1) - (d+2)*S2(d+2));
term4  = (12 * gamma(d+1) * pi_d2) / (3^(d/2) * gamma_d2) * ((a / sqrt3)*S1(d) - S1(d-1));
term5  = (24 * gamma(d+1) * pi_d2) / (3^(d/2) * gamma_d2) * (d*S2(d) - (a*(d+1)/sqrt3)*S2(d+1));
term6  = (36 * gamma(d+2) * pi_d2) / (3^((d+2)/2)) * ((a*(d+2)/sqrt3)*S3(d+2) - d*S3(d+1));
term7  = (12 * gamma(d+2) * pi_d2) / (3^((d+2)/2) * gamma_d2) * ((a*(d+2)/sqrt3)*S1(d+2) - d*S1(d+1));
term8  = ((5*d + 11)/5) * (6 * pi_d2 * gamma(d+2)) / (3^((d+2)/2) * gamma_d2) * ...
         (S1(d+1) - 2*(d+2)*S2(d+2) + (d+2)*(d+3)*S3(d+3));

% Final result
terms = [term1; term2; term3; term4; term5; term6; term7; term8];
R0 = term1 + term2;
R  = R0 + term3 + term4 + term5 + term6 + term7 + term8;
end
