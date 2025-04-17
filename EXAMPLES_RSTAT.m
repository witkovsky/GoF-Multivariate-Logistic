%% EXAMPLES_RSTAT
%  Estimated PDFs for RSTATs (standardized vs. non-standardized)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '16-Apr-2025 17:20:41'

%% EXAMPLE 2
%  Simulated values of the test statistic Rstat for dimension d = 2

d     = 2;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d2_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d2_a2_n100_Std = Rquantile;
MEAN__d2_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,~,~,~,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d2_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d2_a2_n100_nonStd = Rquantile;
MEAN__d2_a2_n100_nonStd = Rmean;
RINF__d2_a2_nInf_zero   = Rinf_zero(a,d);

% Save ALL current variables
save RSTAT_n100

% Create and save Figure RSTAT_d2_a2_n100
 
data1 = RSTAT_d2_a2_n100_Std;
data2 = RSTAT_d2_a2_n100_nonStd;

figure

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'DisplayName','Rstat d=2 a=2 n=100 Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName','Rstat d=2 a=2 n=100 non-Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title('Rstat with and without standardization, d = 2, a = 2, n = 100');

box(axes1,'on');
hold(axes1,'off');
legend(axes1,'show');

savefig('RSTAT_d2_a2_n100.fig')
close(gcf)

%% EXAMPLE 3
%  Simulated values of the test statistic Rstat for dimension d = 3

d     = 3;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d3_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d3_a2_n100_Std = Rquantile;
MEAN__d3_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,~,~,~,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d3_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d3_a2_n100_nonStd = Rquantile;
MEAN__d3_a2_n100_nonStd = Rmean;
RINF__d3_a2_nInf_zero   = Rinf_zero(a,d);

% Save ALL current variables
save RSTAT_n100

% Create and save Figure RSTAT_d3_a2_n100
 
data1 = RSTAT_d3_a2_n100_Std;
data2 = RSTAT_d3_a2_n100_nonStd;

figure

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'DisplayName','Rstat d=3 a=2 n=100 Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName','Rstat d=3 a=2 n=100 non-Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title('Rstat with and without standardization, d = 3, a = 2, n = 100');

box(axes1,'on');
hold(axes1,'off');
legend(axes1,'show');

savefig('RSTAT_d3_a2_n100.fig')
close(gcf)

%% EXAMPLE 4
%  Simulated values of the test statistic Rstat for dimension d = 4

d     = 4;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d4_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d4_a2_n100_Std = Rquantile;
MEAN__d4_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,~,~,~,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d4_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d4_a2_n100_nonStd = Rquantile;
MEAN__d4_a2_n100_nonStd = Rmean;
RINF__d4_a2_nInf_zero   = Rinf_zero(a,d);

% Save ALL current variables
save RSTAT_n100

% Create and save Figure RSTAT_d4_a2_n100
 
data1 = RSTAT_d4_a2_n100_Std;
data2 = RSTAT_d4_a2_n100_nonStd;

figure

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'DisplayName','Rstat d=4 a=2 n=100 Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName','Rstat d=4 a=2 n=100 non-Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title('Rstat with and without standardization, d = 4, a = 2, n = 100');

box(axes1,'on');
hold(axes1,'off');
legend(axes1,'show');

savefig('RSTAT_d4_a2_n100.fig')
close(gcf)

%% EXAMPLE 5
%  Simulated values of the test statistic Rstat for dimension d = 5

d     = 5;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d5_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d5_a2_n100_Std = Rquantile;
MEAN__d5_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,~,~,~,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d5_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d5_a2_n100_nonStd = Rquantile;
MEAN__d5_a2_n100_nonStd = Rmean;
RINF__d5_a2_nInf_zero   = Rinf_zero(a,d);

% Save ALL current variables
save RSTAT_n100

% Create and save Figure RSTAT_d5_a2_n100
 
data1 = RSTAT_d5_a2_n100_Std;
data2 = RSTAT_d5_a2_n100_nonStd;

figure

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'DisplayName','Rstat d=5 a=2 n=100 Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName','Rstat d=5 a=2 n=100 non-Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title('Rstat with and without standardization, d = 5, a = 2, n = 100');

box(axes1,'on');
hold(axes1,'off');
legend(axes1,'show');

savefig('RSTAT_d5_a2_n100.fig')
close(gcf)

%% EXAMPLE 6
%  Simulated values of the test statistic Rstat for dimension d = 6

d     = 6;
N     = 10000;
a     = 2;
alpha = 0.05;
n     = 100;
kMax  = 100;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d6_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d6_a2_n100_Std = Rquantile;
MEAN__d6_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d6_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d6_a2_n100_nonStd = Rquantile;
MEAN__d6_a2_n100_nonStd = Rmean;
RINF__d6_a2_nInf_zero   = Rinf_zero(a,d);

% Save ALL current variables
save RSTAT_n100

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d6_a2_n100_Std;
data2 = RSTAT_d6_a2_n100_nonStd;

figure

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'DisplayName','Rstat d=6 a=2 n=100 Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName','Rstat d=6 a=2 n=100 non-Standardized',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title('Rstat with and without standardization, d = 6, a = 2, n = 100');

box(axes1,'on');
hold(axes1,'off');
legend(axes1,'show');

savefig('RSTAT_d6_a2_n100.fig')
close(gcf)
