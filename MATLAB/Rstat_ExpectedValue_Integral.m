function R = Rstat_ExpectedValue_Integral(d, a, upperLimit)
% RSTAT_EXPECTEDVALUE_INTEGRAL Computes the expected value E[R_{inf,a}] by
%   numerical integration (see proof of Theorem 8 in Popović, Mijanović, 
%   and Witkovský (2025).). This function is based on the integral
%   representation derived from the asymptotic behavior of the multivariate
%   logistic characteristic function.
%
% SYNTAX:
%   R = Rstat_ExpectedValue_Integral(d, a, upperLimit)
%
% INPUT:
%   d - dimension of the multivariate logistic distribution (positive integer)
%   a - positive scalar weight parameter
%   upperLimit - reasonably large upper limit of integration to capture
%       tail decay of exp(-a * r) 
%
% OUTPUT:
%   R - value of E[R_{inf,a}] computed via numerical integration
%
% INTEGRAL FORMULA:
%   Let f(r) = (sqrt(3) * r) / sinh(sqrt(3) * r)
%   Then E[R_{inf,a}]  is given by:
%
%     E[R_{inf,a}]  = C_d * ∫₀^∞ G(r) * r^{d-1} * exp(-a * r) dr
%     where C_d = 2 * π^{d/2} / Γ(d/2)
%
%   and G(r) is defined as:
%     G(r) = 1 - f(r)^2 + r f'(r) (f''(r) + f(r)) + r f(r) f'(r) + ((5d+11)/5) * (r^2 / 4) * f'(r)^2
%
% DEPENDENCIES:
%   Standard MATLAB functions only
%
% EXAMPLES:
%   Rinf = Rstat_ExpectedValue_Integral(2, 2)
%   Rinf = Rstat_ExpectedValue_Integral(2, 2, 20)
%   Rinf = Rstat_ExpectedValue_Integral(3, 2)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '26-Apr-2025 11:04:18'

%% Input validation
if nargin < 3 || isempty(upperLimit), upperLimit = 100; end
if nargin < 2 || isempty(a), a = 2; end
if nargin < 1 || isempty(d), d = 2; end
if d <= 0 || mod(d,1) ~= 0
    error('Dimension d must be a positive integer.');
end
if a <= 0
    error('Parameter a must be positive.');
end

%% Define auxiliary functions
sqrt3 = sqrt(3);
f_r = @(r) sqrt3 .* r ./ sinh(sqrt3 .* r);

f_prime = @(r) sqrt3 ./ sinh(sqrt3 .* r) - (3 .* r .* cosh(sqrt3 .* r)) ./ sinh(sqrt3 .* r).^2;

f_double_prime = @(r) ...
    (6 .* sqrt3 .* r .* cosh(sqrt3 .* r).^2) ./ sinh(sqrt3 .* r).^3 ...
    - (3 .* sqrt3 .* r) ./ sinh(sqrt3 .* r) ...
    - (6 .* cosh(sqrt3 .* r)) ./ sinh(sqrt3 .* r).^2;

%% Define safe integrand to avoid NaNs for large r
integrand = @(r) arrayfun(@(x) ...
    safe_integrand(x, f_r, f_prime, f_double_prime, d, a), r);

%% Coefficient from hyperspherical integration
coeff = 2 * pi^(d / 2) / gamma(d / 2);

%% Compute the integral using adaptive quadrature with upper bound
  
% upperLimit is reasonably large to capture tail decay of exp(-a * r)
R = coeff * integral(integrand, 1e-10, upperLimit, 'RelTol', 1e-10, 'AbsTol', 1e-12);
end

%% Subfunction: stable evaluation of integrand
function val = safe_integrand(r, f_r, f_prime, f_double_prime, d, a)
    try
        fr = f_r(r);
        fp = f_prime(r);
        fpp = f_double_prime(r);

        G = 1 - fr^2 + r * fp * (fpp + fr) + r * fr * fp + ((5 * d + 11) / 5) * (r^2 / 4) * fp^2;
        val = G * r^(d - 1) * exp(-a * r);

        if isnan(val) || isinf(val)
            val = 0;
        end
    catch
        val = 0;
    end
end