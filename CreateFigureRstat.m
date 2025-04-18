function CreateFigureRstat(data1, data2, d, a, n)
%CREATEFIGURE(data1, data2)
%  DATA1:  histogram data
%  DATA2:  histogram data

%  Auto-generated by MATLAB on 18-Apr-2025 09:53:06

% Create figure
figure

% Create axes
axes1 = axes;
hold(axes1,'on');

hold(axes1,'on');
% Create histogram
histogram(data1,'DisplayName',['Rstat d=',num2str(d),' a=',num2str(a),' n=',num2str(n),' Standardized'],...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create histogram
histogram(data2,'DisplayName',['Rstat d=',num2str(d),' a=',num2str(a),' n=',num2str(n),' non-Standardized'],...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('PDF');

% Create xlabel
xlabel('Rstat');

% Create title
title(['Rstat with and without standardization, d = ',num2str(d),' a = ',num2str(a),' n = ',num2str(n)]);

box(axes1,'on');
hold(axes1,'off');
% Create legend
legend(axes1,'show');