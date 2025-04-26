%% RSTAT_CONVERGENCE_SAMPLE-SIZE
%  Estimated PDFs for RSTATs (non-standardized) for different n

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '17-Apr-2025 17:20:41'

clear
close all

%% Set the parameters for simulations
N     = 1000;
d     = 2;
a     = 2;
kMax  = 100;
alpha = 0.05;
n     = [10, 20, 50, 100, 200, 500 1000];
nsz   = length(n);

%% Simulated values of the NON-STANDARDIZED statistic Rstat for ALL chosen sample sizes n

% Non-Standardized Rstat
isStandardized = false;
if isStandardized
    STDname = 'STD';
else
    STDname = 'nonSTD';
end

DATAname  = ['DATA_d',num2str(d),'_a',num2str(a),'_',STDname];
eval([DATAname,' = zeros(',num2str(N),',',num2str(nsz),');'])
for i = 1:nsz
    nid = n(i);
    [Rquantile,Rmean,Rstats] = RCriticalValues(d,N,a,alpha,nid,kMax,isStandardized);

    RSTATname = ['RSTAT_d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];
    QSTATname = ['QSTAT_d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];
    MEANname  = ['MEAN__d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];

    eval([RSTATname,' = reshape(Rstats(1,1,:),N,1);'])
    eval([QSTATname,' = Rquantile;'])
    eval([MEANname,' = Rmean;'])
    eval([DATAname,'(:,',num2str(i),') = ',RSTATname])
end

% Plot figure of estimated PDFs (histograms) of RSTATs for ALL sample sizes n
CreateFigureRstatConv(eval(DATAname), d, a, n, isStandardized)

% Save Figure
FIGname  = ['FIG_RSTAT_CONV_d',num2str(d),'_a',num2str(a),'_',STDname,'.fig'];
savefig(FIGname)

% Save the current variables
FILEname  = ['RSTAT_CONVERGENCE_DATA_d',num2str(d),'_a',num2str(a),'_',STDname];
save(FILEname)

%% Simulated values of the STANDARDIZED statistic Rstat for ALL chosen sample sizes n

% Standardized Rstat
isStandardized = true;
if isStandardized
    STDname = 'STD';
else
    STDname = 'nonSTD';
end

DATAname  = ['DATA_d',num2str(d),'_a',num2str(a),'_',STDname];
eval([DATAname,' = zeros(',num2str(N),',',num2str(nsz),');'])
for i = 1:nsz
    nid = n(i);
    [Rquantile,Rmean,Rstats] = RCriticalValues(d,N,a,alpha,nid,kMax,isStandardized);

    RSTATname = ['RSTAT_d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];
    QSTATname = ['QSTAT_d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];
    MEANname  = ['MEAN__d',num2str(d),'_a',num2str(a),'_n',num2str(nid),'_',STDname];

    eval([RSTATname,' = reshape(Rstats(1,1,:),N,1);'])
    eval([QSTATname,' = Rquantile;'])
    eval([MEANname,' = Rmean;'])
    eval([DATAname,'(:,',num2str(i),') = ',RSTATname])
end

% Plot figure of estimated PDFs (histograms) of RSTATs for ALL sample sizes n
CreateFigureRstatConv(eval(DATAname), d, a, n, isStandardized)

% Save Figure
FIGname  = ['FIG_RSTAT_CONV_d',num2str(d),'_a',num2str(a),'_',STDname,'.fig'];
savefig(FIGname)

% Save the current variables
FILEname  = ['RSTAT_CONVERGENCE_DATA_d',num2str(d),'_a',num2str(a),'_',STDname];
save(FILEname)
