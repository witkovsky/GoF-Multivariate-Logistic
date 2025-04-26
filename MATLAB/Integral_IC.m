function IC = Integral_IC(m, a, method, upperLimit)
% INTEGRAL_IC Computes the integral I_C(m, a) based on Lemma 9 in Popović,
% Mijanović and Witkovský (2025).  
%
% The integral is given by:
%   I_C(m, a) = ∫₀^∞ [r^m * coth^2(sqrt(3)*r) / sinh^2(sqrt(3)*r)] * exp(-a*r) dr
%
% It can be computed:
%   - By direct numerical integration
%   - By an expression involving the previously computed I_A and I_B integrals
%
% SYNTAX:
%   IC = Integral_IC(m, a)
%   IC = Integral_IC(m, a, method, upperLimit)
%
% INPUT:
%   m      - integer >= 4 (degree)
%   a      - positive real scalar (weight parameter)
%   method - (optional) 'numerical' (default) or 'zeta'
%            'numerical' - direct numerical integration
%            'zeta'      - uses computation based on I_A and I_B
%   upperLimit  - upper limit of integration of J(mu,beta), default value is
%            upperLimit = 100. 
%
% OUTPUT:
%   IC - computed value of I_C(m, a)
%
% DEPENDENCIES:
%   Integral_IA.m, Integral_IB.m
%
% EXAMPLES:
%   Integral_IC(5, 2)
%   Integral_IC(5, 2, 'zeta')
%   Integral_IC(5, 2, 'integral')
%   Integral_IC(5, 2, 'integral', 20)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver. '26-Apr-2025 11:48:16'


%% Input validation
if nargin < 3 || isempty(method), method = 'zeta'; end
if nargin < 4 || isempty(upperLimit), upperLimit = 100; end

if nargin < 4
    method = 'integral';
end

if a <= 0
    error('Parameter a must be positive.');
end

sqrt3 = sqrt(3);

switch lower(method)
    case 'integral'
        % Direct numerical integral integration
        integrand = @(r) safe_integrand(r, m, a, sqrt3);
        IC = integral(integrand, 0, upperLimit, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    case 'zeta'
        % Using Lemma 9 expression based on I_A and I_B
        const1 = 2 / 3^((m+3)/2);
        const2 = 1 / (3*sqrt3);
        term1  = a * Integral_J(m+1,a/sqrt3,'zeta') - m * Integral_J(m,a/sqrt3,'zeta');
        term2  = m * Integral_IB(m-1,a,'zeta') - a * Integral_IB(m,a,'zeta');
        term3  = Integral_IA(m,a, 'zeta');
        IC = const1 *  term1 + const2 * term2 + term3;

    otherwise
        error('Unknown method. Choose either "integral" or "zeta".');
end
end

%% Subfunction for safe evaluation of integrand
function val = safe_integrand(x, m, a, sqrt3)
if x <= eps
    x = eps;
    val = exp(m*log(x) + 2*log(coth(sqrt3*x)) ...
        - 2*log(sinh(sqrt3*x)) - a * x);
else
    val = (x.^m .* coth(sqrt3*x).^2 ./ sinh(sqrt3*x).^2) .* exp(-a*x);
end
end