function result = SumS3(m, a, N, M)
% SUMS3_EM Efficiently evaluates the infinite sum:
%   S3(m, a) = sum_{n=1}^{∞} n^2 * [2^{-m} * zeta(m+1, n + a/(2√3)) - (2n + a/√3)^{-(m+1)}]
%
% using the Euler–Maclaurin summation formula for fast convergence.
%
% SYNTAX:
%   result = SumS3(m, a)
%   result = SumS3(m, a, N)
%   result = SumS3(m, a, N, M)
%
% INPUT:
%   m - non-negative integer exponent
%   a - positive real parameter
%   N - number of explicit terms in sum (default: 100)
%   M - number of Bernoulli correction terms (default: 2)
%
% OUTPUT:
%   result - approximation of the infinite sum
%
% EXAMPLES:
%   SumS3(2, 1)
%   SumS3(4, 0.5, 100)
%

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-May-2024 16:29:07

if nargin < 3, N = 1000; end
if nargin < 4, M = 2; end

if M > 2
    warning('M - number of Bernoulli correction terms is too large. Set M = 2 !')
    M = 2;
end

% Bernoulli numbers B_{2k}, k = 1:8
B = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510];

sqrt3 = sqrt(3);
f = @(x) x.^2 .* (2^(-m) .* HurwitzZeta(m+1, x + a/(2*sqrt3)) - (2*x + a/sqrt3).^(-(m+1)));

% Step 1: Direct summation for first N-1 terms
n = 1:N-1;
sum_main = sum(f(n));

% Step 2: Integral approximation of tail
if M > 0
    x = N;
    F_tail = integral(f, x, Inf, 'AbsTol', 1e-12);
    F_corr = -0.5 * f(x);

    % Step 3: Euler–Maclaurin correction with numerical derivatives
    deriv_sum = 0;
    for k = 1:M
        B2k = B(k);
        coeff = B2k / factorial(2*k);
        deriv = nth_derivative(f, x, 2*k - 1);
        deriv_sum = deriv_sum + coeff * deriv;
    end
    result = sum_main + F_tail + F_corr + deriv_sum;
else
    result = sum_main;
end
end

%% Finite difference derivative
function d = nth_derivative(f, x, n)
h = 1e-4;
c = central_diff_coeffs(n);
m = floor(length(c)/2);
d = 0;
for i = 1:length(c)
    d = d + c(i) * f(x + (i - m - 1)*h);
end
d = d / h^n;
end

%% Central difference coefficients (up to n = 15)
function c = central_diff_coeffs(n)
switch n
    case 1,  c = [1/2, 0, -1/2];
    case 2,  c = [-1, 2, -1];
    case 3,  c = [-1/2, 1, 0, -1, 1/2];
    case 4,  c = [1, -4, 6, -4, 1];
    case 5,  c = [1, -5, 10, 0, -10, 5, -1]/2;
    case 6,  c = [-1, 6, -15, 20, -15, 6, -1];
    case 7,  c = [-1, 7, -21, 35, 0, -35, 21, -7, 1]/2;
    case 8,  c = [1, -8, 28, -56, 70, -56, 28, -8, 1];
    case 9,  c = [1, -9, 36, -84, 126, 0, -126, 84, -36, 9, -1]/2;
    case 10, c = [-1, 10, -45, 120, -210, 252, -210, 120, -45, 10, -1];
    case 11, c = [-1, 11, -55, 165, -330, 462, 0, -462, 330, -165, 55, -11, 1]/2;
    case 12, c = [1, -12, 66, -220, 495, -792, 924, -792, 495, -220, 66, -12, 1];
    case 13, c = [1, -13, 78, -286, 715, -1287, 1716, 0, -1716, 1287, -715, 286, -78, 13, -1]/2;
    case 14, c = [-1, 14, -91, 364, -1001, 2002, -3003, 3432, -3003, 2002, -1001, 364, -91, 14, -1];
    case 15, c = [-1, 15, -105, 455, -1365, 3003, -5005, 6435, 0, -6435, 5005, -3003, 1365, -455, 105, -15, 1]/2;
    otherwise
        error('Finite difference coefficients not implemented for n > 15');
end
end