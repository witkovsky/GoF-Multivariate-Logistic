function IA = Integral_IA(m, a, method, upperLimit)
% INTEGRAL_IA Computes the integral I_A(m, a) based on Lemma 7 in Popović,
% Mijanović and Witkovský (2025).  
%
% The integral is given by:
%   I_A(m, a) = ∫₀^∞ [r^m / sinh^2(sqrt(3) * r)] * exp(-a * r) dr
%
% It can be computed:
%   - By direct numerical integration
%   - By an expression involving J(mu, beta) integrals from Lemma 3
%
% SYNTAX:
%   IA = Integral_IA(m, a)
%   IA = Integral_IA(m, a, method, upperLimit)
%
% INPUT:
%   m      - integer >= 2 (degree)
%   a      - positive real scalar (weight parameter)
%   method - (optional) 'integral' or 'zeta' (default) 
%            'integral' - direct numerical integration
%            'zeta'      - uses Integral_J for semi-analytical computation
%   upperLimit  - upper limit of integration of J(mu,beta), default value is
%            upperLimit = 100. 
%
% OUTPUT:
%   IA - computed value of I_A(m, a)
%
% DEPENDENCIES:
%   Integral_J.m (for 'zeta' method)
%
% EXAMPLES:
%   Integral_IA(3, 1)
%   Integral_IA(3, 1, 'zeta')
%   Integral_IA(3, 1, 'integral')
%   Integral_IA(3, 1, 'integral', 20)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver. '26-Apr-2025 11:48:16'

%% Input validation
if nargin < 3 || isempty(method), method = 'zeta'; end
if nargin < 4 || isempty(upperLimit), upperLimit = 100; end

if m < 2 
    method = 'integral';
end

if a <= 0
    error('Parameter a must be positive.');
end

sqrt3 = sqrt(3);

switch lower(method)

    case 'integral'
        % Direct numerical integration
        integrand = @(r) safe_integrand(r, m, a, sqrt3);
        IA = integral(integrand, 0, upperLimit, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    case 'zeta'
        % Using Lemma 7 formula
        beta = a / sqrt3;
        term1 = m * Integral_J(m, beta, 'zeta');
        term2 = beta * Integral_J(m + 1, beta, 'zeta');
        IA = 3^(-(m+1)/2) * (term1 - term2);

    otherwise
        error('Unknown method. Choose either "numerical" or "zeta".');
end
end

%% Subfunction for safe evaluation of integrand
function val = safe_integrand(x, m, a, sqrt3)
if x <= eps
    x = eps;
    val = exp( m.*log(x) - 2 * log(sinh(sqrt3 * x)) -a.* x );
else
    val = (x.^m ./ sinh(sqrt3 * x).^2) .* exp(-a * x);
end
end
