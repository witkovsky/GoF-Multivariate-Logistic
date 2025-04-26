function J = Integral_J(m, a, method, upperLimit)
% INTEGRAL_J Computes the integral J(m, a) defined as:
%   J(m, a) = ∫₀^∞ x^{m-1} exp(-a*x) coth(x) dx
% according to Lemma 3 in Popović, Mijanović and Witkovský (2025).
%
% SYNTAX:
%   J = Integral_J(m, a)
%   J = Integral_J(m, a, method)
%
% INPUT:
%   m      - real scalar, Re(m) > 1
%   a      - real scalar, a > 0
%   method - (optional) 'zeta' (default) or 'integral'
%            'zeta'      - uses closed-form formula with Hurwitz zeta function
%            'integral' - computes directly via numerical integration
%   upperLimit  - upper limit of integration of J(m,a), default value is
%            upperLimit = 100. 
%
% OUTPUT:
%   J - value of the integral
%
% DEPENDENCIES:
%   HurwitzZeta.m (if using 'zeta' method)
%
% EXAMPLES:
%   Integral_J(2, 1)
%   Integral_J(2, 1, 'zeta')
%   Integral_J(2, 1, 'integral')
%   Integral_J(2, 1, 'integral', 20)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver. '26-Apr-2025 11:48:16'

%% Input checks
if nargin < 3 || isempty(method), method = 'zeta'; end
if nargin < 4 || isempty(upperLimit), upperLimit = 100; end

if m <= 1
    method = 'integral';
end

if a <= 0
    error('Parameter a must be positive.');
end

switch lower(method)
    
    case 'zeta'
        % Method based on Lemma 3 formula
        z = a / 2;
        J = gamma(m) .* (2.^(1 - m) .* HurwitzZeta(m, z) - a.^(-m));

    case 'integral'
        % Numerical integration method
        integrand = @(x) safe_integrand(x, m, a);
        J = integral(integrand, 0, upperLimit, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    otherwise
        error('Unknown method. Choose "zeta" or "numerical".');
end
end

%% Subfunction for safe evaluation of the integrand
function val = safe_integrand(x, m, a)
if x == 0
    x = eps;
    val = exp( (m-1).*log(x) + log(coth(x)) - a .* x);
else
    val = x.^(m-1) .* coth(x) .* exp(- a .* x);
end
end