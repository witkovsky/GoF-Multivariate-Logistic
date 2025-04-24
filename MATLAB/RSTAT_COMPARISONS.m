%% RSTAT_COMPARISONS_STANDARDIZATION
%  Estimated PDFs for RSTATs (standardized vs. non-standardized)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '17-Apr-2025 17:20:41'

clear
close all

%% Set the parameters for simulations
N     = 10000;
a     = 2;
n     = 100;

kMax  = 100;
alpha = 0.05;

%% EXAMPLE 2
%  Simulated values of the test statistic Rstat for dimension d = 2

% dimension
d     = 2;

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

% Create and save Figure RSTAT_d2_a2_n100
 
data1 = RSTAT_d2_a2_n100_Std;
data2 = RSTAT_d2_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 3
%  Simulated values of the test statistic Rstat for dimension d = 3
clear
close all

% dimension
d     = 3;

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

% Create and save Figure RSTAT_d3_a2_n100
 
data1 = RSTAT_d3_a2_n100_Std;
data2 = RSTAT_d3_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 4
%  Simulated values of the test statistic Rstat for dimension d = 4
clear
close all

% dimension
d     = 4;

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

% Create and save Figure RSTAT_d4_a2_n100
 
data1 = RSTAT_d4_a2_n100_Std;
data2 = RSTAT_d4_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 5
%  Simulated values of the test statistic Rstat for dimension d = 5
clear
close all

% dimension
d     = 5;

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

% Create and save Figure RSTAT_d5_a2_n100
 
data1 = RSTAT_d5_a2_n100_Std;
data2 = RSTAT_d5_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 6
%  Simulated values of the test statistic Rstat for dimension d = 6
clear
close all

% dimension
d     = 6;

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

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d6_a2_n100_Std;
data2 = RSTAT_d6_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 7
%  Simulated values of the test statistic Rstat for dimension d = 7
clear
close all

% dimension
d     = 7;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d7_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d7_a2_n100_Std = Rquantile;
MEAN__d7_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d7_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d7_a2_n100_nonStd = Rquantile;
MEAN__d7_a2_n100_nonStd = Rmean;
RINF__d7_a2_nInf_zero   = Rinf_zero(a,d);

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d7_a2_n100_Std;
data2 = RSTAT_d7_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)


%% EXAMPLE 8
%  Simulated values of the test statistic Rstat for dimension d = 8
clear
close all

% dimension
d     = 8;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d8_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d8_a2_n100_Std = Rquantile;
MEAN__d8_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d8_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d8_a2_n100_nonStd = Rquantile;
MEAN__d8_a2_n100_nonStd = Rmean;
RINF__d8_a2_nInf_zero   = Rinf_zero(a,d);

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d8_a2_n100_Std;
data2 = RSTAT_d8_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)

%% EXAMPLE 9
%  Simulated values of the test statistic Rstat for dimension d = 10
clear
close all

% dimension
d     = 9;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d9_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d9_a2_n100_Std = Rquantile;
MEAN__d9_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d9_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d9_a2_n100_nonStd = Rquantile;
MEAN__d9_a2_n100_nonStd = Rmean;
RINF__d9_a2_nInf_zero   = Rinf_zero(a,d);

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d9_a2_n100_Std;
data2 = RSTAT_d9_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)

% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)


%% EXAMPLE 10
%  Simulated values of the test statistic Rstat for dimension d = 10
clear
close all

% dimension
d     = 10;

% Rsrat calculated from standardized values Z = \ha{Sigma}^{-1/2}(X-\ha{mu})

isStandrdized = true;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d10_a2_n100_Std = reshape(Rstats(1,1,:),N,1);
QUANT_d10_a2_n100_Std = Rquantile;
MEAN__d10_a2_n100_Std = Rmean;

% Rsrat calculated from non-standardized values X ~ ML(0,I)

isStandrdized = false;
[Rquantile,Rmean,d,N,a,alpha,n,kMax,Rstats] = RCriticalValues(d,N,a,alpha,n,kMax,isStandrdized);
disp(Rquantile)

RSTAT_d10_a2_n100_nonStd = reshape(Rstats(1,1,:),N,1);
QUANT_d10_a2_n100_nonStd = Rquantile;
MEAN__d10_a2_n100_nonStd = Rmean;
RINF__d10_a2_nInf_zero   = Rinf_zero(a,d);

% Create and save Figure RSTAT_d6_a2_n100
 
data1 = RSTAT_d10_a2_n100_Std;
data2 = RSTAT_d10_a2_n100_nonStd;
CreateFigureRstat(data1, data2, d, a, n)

% Save Figure
FIGname  = ['FIG_RSTAT_COMP_d',num2str(d),'_a',num2str(a),'_n',num2str(n),'.fig'];
savefig(FIGname)

% close(gcf)
% Save the current variables
FILEname  = ['RSTAT_COMPARISON_DATA_d',num2str(d),'_a',num2str(a),'_n',num2str(n)];
save(FILEname)