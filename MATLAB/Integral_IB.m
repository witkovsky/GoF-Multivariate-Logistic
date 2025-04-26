function IB = Integral_IB(m, a, method, upperLimit)
% INTEGRAL_IB Computes the integral I_B(m, a) based on Lemma 8 in Popović,
% Mijanović and Witkovský (2025).  
%
% The integral is given by:
%   I_B(m, a) = ∫₀^∞ [r^m * coth(sqrt(3)*r) / sinh^2(sqrt(3)*r)] * exp(-a*r) dr
%
% It can be computed:
%   - By direct numerical integration (first part of Lemma 8)
%   - By an expression involving the previously computed I_A integrals
%
% SYNTAX:
%   IB = Integral_IB(m, a)
%   IB = Integral_IB(m, a, method, upperLimit)
%
% INPUT:
%   m      - integer >= 3 (degree)
%   a      - positive real scalar (weight parameter)
%   method - (optional) 'integral' or 'zeta' (default)
%            'integral' - direct numerical integration of the definition
%            'zeta'     - uses semi-analytical computation based on I_A
%   upperLimit  - upper limit of integration of J(mu,beta), default value is
%            upperLimit = 100. 
%
% OUTPUT:
%   IB - computed value of I_B(m, a)
%
% DEPENDENCIES:
%   Integral_IA.m (for 'zeta' method)
%
% EXAMPLES:
%   Integral_IB(4, 2)
%   Integral_IB(4, 2, 'zeta')
%   Integral_IB(4, 2, 'integral')
%   Integral_IB(4, 2, 'integral', 20)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver. '26-Apr-2025 11:48:16'

%% Input validation
if nargin < 3 || isempty(method), method = 'zeta'; end
if nargin < 4 || isempty(upperLimit), upperLimit = 100; end

if m < 3 
    method = 'integral';
end

if a <= 0
    error('Parameter a must be positive.');
end

sqrt3 = sqrt(3);

switch lower(method)
    case 'integral'
        % Direct numerical integration of the definition
        integrand = @(r) safe_integrand(r, m, a, sqrt3);
        IB = integral(integrand, 0, upperLimit, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    case 'zeta'
        % Using Lemma 8 expression based on I_A
        IA_m_minus_1 = Integral_IA(m-1, a, 'zeta');
        IA_m         = Integral_IA(m, a, 'zeta');
        IB = (m / (2 * sqrt3)) * IA_m_minus_1 - (a / (2 * sqrt3)) * IA_m;

    otherwise
        error('Unknown method. Choose either "numerical" or "zeta".');
end
end

%% Subfunction for safe evaluation of integrand
function val = safe_integrand(x, m, a, sqrt3)
if x <= eps
    x = eps;
    val = exp(m * log(x) + log(coth(sqrt3 * x)) ...
        - 2 * log(sinh(sqrt3 * x)) - a * x);
else
    val = (x.^m .* coth(sqrt3 * x) ./ sinh(sqrt3 * x).^2) .* exp(-a * x);
end
end