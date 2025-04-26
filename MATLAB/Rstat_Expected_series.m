function Rinf = Rstat_Expected_series(d, a, method, upperLimit)
% RSTAT_EXPECTED_SERIES Computes the expected value E[R_{inf,a}] based on
%   Theorem 8 in Popović, Mijanović, and Witkovský (2025). The computation
%   combines the integrals J, I_A, I_B, and I_C according to their
%   definitions in Lemma 3, and Lemmas 7-9.
%
% SYNTAX:
%   Rinf = Rstat_Expected_series(d, a)
%   Rinf = Rstat_Expected_series(d, a, method, upperLimit)
%
% INPUT:
%   d      - integer, dimension (d >= 3)
%   a      - positive real scalar, weight parameter
%   method - (optional) 'zeta' (default) or 'integral'
%            'zeta'      - uses semi-analytical expressions
%            'integral'  - computes by direct numerical integration
%   upperLimit - reasonably large upper limit of integration to capture
%       tail decay of exp(-a * r) 
%
% OUTPUT:
%   Rinf - computed expected value of R_inf,a, E[R_{inf,a}]
%
% DEPENDENCIES:
%   Integral_J.m, Integral_IA.m, Integral_IB.m, Integral_IC.m
%
% EXAMPLES:
%   Rinf = Rstat_Expected_series(3, 2)
%   Rinf = Rstat_Expected_series(3, 2, 'zeta')
%   Rinf = Rstat_Expected_series(3, 2, 'integral')
%   Rinf = Rstat_Expected_series(3, 2, 'integral', 20)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '26-Apr-2025 11:04:18'

%% Input validation
if nargin < 3 || isempty(method), method = 'zeta'; end
if nargin < 4 || isempty(upperLimit), upperLimit = 100; end

if nargin < 4
    method = 'integral';
end

if a <= 0
    error('Parameter a must be positive.');
end

%% Constants
sqrt3 = sqrt(3);
pi_d2 = pi^(d/2);
gamma_d2 = gamma(d/2);

%% Compute required integrals
IA_d_plus1 = Integral_IA(d+1, a, method, upperLimit);
IA_d       = Integral_IA(d, a, method, upperLimit);
IA_d_minus1 = Integral_IA(d-1, a, method, upperLimit);

IB_d_plus2 = Integral_IB(d+2, a, method, upperLimit);
IB_d_plus1 = Integral_IB(d+1, a, method, upperLimit);
IB_d       = Integral_IB(d, a, method, upperLimit);

IC_d_plus3 = Integral_IC(d+3, a, method, upperLimit);
IC_d_plus2 = Integral_IC(d+2, a, method, upperLimit);
IC_d_plus1 = Integral_IC(d+1, a, method, upperLimit);

%% Assemble the expression according to Theorem 8

prefactor = (2 * pi_d2 * gamma(d)) / (gamma_d2 * a^d);

sum_terms = (15*d + 93)/10 * IA_d_plus1 ...
           + 3*a * IA_d ...
           - 3*d * IA_d_minus1 ...
           - (5*d + 31)/10 * 6*sqrt3 * IB_d_plus2 ...
           + 6*d*sqrt3 * IB_d ...
           - 6*a*sqrt3 * IB_d_plus1 ...
           + 9*a * IC_d_plus2 ...
           - 9*d * IC_d_plus1 ...
           + (5*d + 11)/10 * 9 * IC_d_plus3;

Rinf = prefactor + (pi_d2 / gamma_d2) * sum_terms;

end
