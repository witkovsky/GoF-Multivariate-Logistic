function R = Rinf_zero(a, d)
% RINF_ZERO Computes E[R_inf,a] under the null hypothesis H0
% that the sample comes from a multivariate logistic distribution
% with mean zero and identity covariance matrix.
%
% This matches the simplified representation used in Theorem 6.
%
% SYNTAX:
%   R = Rinf_zero(a, d)
%
% INPUT:
%   a - positive scalar, weight parameter
%   d - positive integer, dimension
%
% OUTPUT:
%   R - expected value E[R_inf,a] under H0
%
% FORMULA:
%   R = C1 - C2 * (zeta(d+1, z) - z * zeta(d+2, z))
%   where:
%       z    = a / (2 * sqrt(3))
%       C1   = 2 * pi^{d/2} * gamma(d) / (a^d * gamma(d/2))
%       C2   = 2^{1-d} * pi^{d/2} * gamma(d+2) / (3^{d/2} * gamma(d/2))
%
% DEPENDENCIES:
%   HurwitzZeta.m  (vectorized implementation of zeta(s, q))
%
% EXAMPLE:
%   Rinf_zero(2, 3)

% Input validation
if nargin < 2
    error('Two input arguments required: a (scalar), d (dimension).');
end
if a <= 0 || d <= 0 || mod(d, 1) ~= 0
    error('Input a must be > 0 and d must be a positive integer.');
end

% Constants
sqrt3 = sqrt(3);
z = a / (2 * sqrt3);

C1 = (2 * pi^(d / 2) * gamma(d)) / (a^d * gamma(d / 2));
C2 = (2^(1 - d) * pi^(d / 2) * gamma(d + 2)) / (3^(d / 2) * gamma(d / 2));

% Hurwitz zeta evaluations
Z1 = HurwitzZeta(d + 1, z);
Z2 = HurwitzZeta(d + 2, z);

% Final expected value under H0
R = C1 - C2 * (Z1 - z * Z2);
end
