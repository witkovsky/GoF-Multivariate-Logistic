function z = HurwitzZeta(s, q)
% HurwitzZeta Numerically evaluates the Hurwitz zeta function \zeta(s, q)
% for complex parameter s with Re(s) > 1 and q > 0.
%
% The Hurwitz zeta function is defined as:
%     \zeta(s, q) = \sum_{n=0}^{\infty} 1 / (q + n)^s
% and generalizes the Riemann zeta function \zeta(s) = \zeta(s, 1).
%
% This implementation uses the Euler–Maclaurin summation formula for
% fast and accurate convergence.
%
% SYNTAX:
%   z = HurwitzZeta(s, q)
%
% INPUT:
%   s - complex scalar with real(s) > 1
%   q - positive real scalar (q > 0)
%
% OUTPUT:
%   z - numerical approximation to \zeta(s, q)
%
% NOTES:
%   - Accurate for a wide range of s and q, especially s not too close to 1
%   - For symbolic evaluations, use MATLAB's hurwitzZeta(s,q) if Symbolic
%     Toolbox is available 
%
% EXAMPLES:
%   HurwitzZeta(2, 1)        % returns pi^2/6
%   HurwitzZeta(3, 1)        % matches Apery's constant
%   HurwitzZeta(2.5, 0.5)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '16-Apr-2025 14:30:43'

%% Parameters
N = 50;            % Number of terms in the initial sum
M = 8;             % Number of Bernoulli correction terms (B_{2k}, k = 1:8)

% Predefined Bernoulli numbers B_{2k}
B = [ 1/6, ...       % B_2
     -1/30, ...      % B_4
      1/42, ...      % B_6
     -1/30, ...      % B_8
      5/66, ...      % B_10
    -691/2730, ...   % B_12
      7/6, ...       % B_14
   -3617/510 ];      % B_16

%% Input checks
if nargin < 2
    error('HurwitzZeta requires two input arguments: s and q');
end
if real(s) <= 1 || q <= 0
    error('HurwitzZeta is defined for real(s) > 1 and q > 0');
end

%% Initial sum: sum_{k=0}^{N-1} (q + k)^(-s)
z = sum((q + (0:N-1)).^(-s));

%% Euler–Maclaurin correction
qN = q + N;
z = z + (qN^(1 - s)) / (s - 1) + 0.5 * qN^(-s);

% Bernoulli correction terms
for k = 1:M
    B2k = B(k);
    term = B2k * gamma(s + 2*k - 1) / factorial(2*k) / qN^(s + 2*k - 1);
    z = z + term;
end
end