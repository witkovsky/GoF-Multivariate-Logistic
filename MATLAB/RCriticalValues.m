function [Rstat_quantile, Rstat_mean, d, N, a, alpha, n, kMax, Rstat_values] = ...
    RCriticalValues(d, N, a, alpha, n, kMax, isStandardized)
% RCRITICALVALUES  Computes critical values of the GOF test statistic R
% by Monte Carlo simulation for given dimensions, weight parameters,
% significance levels, and sample sizes.
%
% SYNTAX:
%   [Rstat_quantile, Rstat_mean, d, N, a, alpha, n, kMax, Rstat_values] = ...
%                 RCriticalValues(d, N, a, alpha, n, kMax, isStandardized)
%
% INPUT:
%   d      - dimension of the multivariate logistic distribution (default: 2)
%   N      - number of Monte Carlo simulations (default: 1000)
%   a      - vector of weight parameters for the test statistic (default: 2)
%   alpha  - vector of significance levels (default: 0.05)
%   n      - vector of sample sizes (default: 100)
%   kMax   - truncation level for the I3-series in Rstat (default: 100)
%   isStandardized - indicator for standardization of the samples (default: true)
%
% OUTPUT:
%   Rstat_quantile - matrix of critical values for each (a, alpha, n)
%   Rstat_mean     - estimated mean of R statistics for each (a, n)
%   Rstat_values   - simulated test statistic values (na x nn x N)
%
% EXAMPLE 1: (Small simulation)
%   d = 2; 
%   N = 100; 
%   a = [1 2];  
%   alpha = [0.05 0.01]; 
%   n = [50 100]; 
%   kMax = 100;
%   [Rstat_quantile, Rstat_mean] = RCriticalValues(d, N, a, alpha, n, kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)
%
% EXAMPLE 2: Bivariate case
%   d     = 2;
%   N     = 1000;
%   a     = 2;
%   alpha = 0.05;
%   n     = 100;
%   kMax  = 100;
%   [Rstat_quantile, Rstat_mean, d, N, a, alpha, n, kMax, Rstat_values] = RCriticalValues(d,N,a,alpha,n,kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)
%
% EXAMPLE 3: Trivariate case
%   d     = 3;
%   N     = 1000;
%   a     = 2;
%   alpha = 0.05;
%   n     = 100;
%   kMax  = 100;
%   [Rstat_quantile,Rstat_mean,d,N,a,alpha,n,kMax,Rstat_values] = RCriticalValues(d,N,a,alpha,n,kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)
%
% EXAMPLE 4: 4-dimensional case
%   d     = 4;
%   N     = 1000;
%   a     = 2;
%   alpha = 0.05;
%   n     = 100;
%   kMax  = 100;
%   [Rstat_quantile,Rstat_mean,d,N,a,alpha,n,kMax,Rstat_values] = RCriticalValues(d,N,a,alpha,n,kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)
% 
% EXAMPLE 5 (Table 1 from the Popović, Mijanović, Witkovský (2025))
%   rng("default")
%   d     = 2;
%   N     = 10000;
%   a     = [0.5 1 1.5 2 2.5 3 3.5 5];
%   alpha = 0.05;
%   n     = [50 100];
%   kMax  = 100;
%   [Rstat_quantile,Rstat_mean,d,N,a,alpha,n,kMax,Rstat_values] = RCriticalValues(d,N,a,alpha,n,kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)
% 
% EXAMPLE 6 (Related to Table 2 from the Popović, Mijanović, Witkovský (2025))
%   rng("default")
%   d     = 3;
%   N     = 10000;
%   a     = [0.5 1 1.5 2 2.5 3 3.5 5];
%   alpha = 0.05;
%   n     = [50 100];
%   kMax  = 100;
%   [Rstat_quantile,Rstat_mean,d,N,a,alpha,n,kMax,Rstat_values] = RCriticalValues(d,N,a,alpha,n,kMax);
%   disp(Rstat_mean) 
%   disp(Rstat_quantile)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.:  '26-Apr-2025 10:02:24'

%% Input handling
narginchk(0, 7);
if nargin < 7 || isempty(isStandardized),  isStandardized = true;  end
if nargin < 6 || isempty(kMax),  kMax = 100;  end
if nargin < 5 || isempty(n),     n = 100;     end
if nargin < 4 || isempty(alpha), alpha = 0.05; end
if nargin < 3 || isempty(a),     a = 2;        end
if nargin < 2 || isempty(N),     N = 1000;     end
if nargin < 1 || isempty(d),     d = 2;        end

%% Preprocessing and initialization
a      = a(:);         % ensure column vector
alpha  = alpha(:);
n      = n(:);
na     = length(a);    % number of weight parameters
nalpha = length(alpha);% number of significance levels
nn     = length(n);    % number of sample sizes

Rstat_quantile = zeros(na, nalpha, nn);   % critical values
Rstat_mean     = zeros(na, nn);           % mean of R stats
Rstat_values   = zeros(na, nn, N);        % full simulation data

mu    = zeros(1, d);       % mean vector
Sigma = eye(d);            % identity covariance matrix

%% Monte Carlo Simulation
for i = 1:na          % loop over weights a(i)
    for j = 1:nn      % loop over sample sizes n(j)
        for kk = 1:N
            % Display simulation progress (can comment out for speed)
            fprintf('Sim %d/%d | d = %d | a = %.3f | n = %d\n', kk, N, d, a(i), n(j));

            % Generate sample from multivariate logistic distribution
            data = randML(n(j), mu, Sigma);

            % Compute test statistic
            Rstat_values(i, j, kk) = Rstat(data, a(i), kMax, isStandardized);
        end

        % Sort and store statistics
        Rstat_sorted = sort(Rstat_values(i, j, :));
        Rstat_values(i, j, :) = Rstat_sorted;

        % Estimate expected value of R-stat
        Rstat_mean(i, j) = mean(Rstat_sorted);

        % Estimate quantiles for all specified significance levels
        for k = 1:nalpha
            Rstat_quantile(i, k, j) = Rstat_sorted(ceil(N * (1 - alpha(k))));
        end
    end
end
end
