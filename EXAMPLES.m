%% EXAMPLES

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.:  '02-May-2024 08:06:49'

%% EXAMPLE 1
%  Simulated values of the test statistic for d = 3

d     = 3;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,TT] = RCriticalValues(d,N,a,alpha,n,kMax);
disp(Rquantile)

TD3 = reshape(TT(1,1,:),N,1);
QD3 = Rquantile
MD3 = Rmean

%% EXAMPLE 1
%  Simulated values of the test statistic for d = 4

d     = 4;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,TT] = RCriticalValues(d,N,a,alpha,n,kMax);
disp(Rquantile)

TD4 = reshape(TT(1,1,:),N,1);
QD4 = Rquantile;
MD4 = Rmean;

%%  REAL EXAMPLE 
%   Simulated values of the test statistic, d = 2, n = 87, a = 2

d     = 2;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 87;
kMax  = 100;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,TT] = RCriticalValues(d,N,a,alpha,n,kMax);
disp(Rquantile)

TD2n87 = reshape(TT(1,1,:),N,1);
QD2n87 = Rquantile;
MD2n87 = Rmean;
