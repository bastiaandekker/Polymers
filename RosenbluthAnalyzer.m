% Script for analyzing the data files produced by rosenbluth.f90
% B.Dekker - ICCP 2014
%%
clear all;

%% Read data from file, separating the data  in a number of data blocks
numDatablocks = 20;
[blockSizeRosen, beadsRosen, r2Rosen] = analyzeRosenbluthFile('RosenbluthCalc.txt', numDatablocks);
[blockSizePerm,  beadsPerm, r2Perm] = analyzeRosenbluthFile('PERMCalc.txt', numDatablocks);


%% Plot
close all;
% Exponential fit
a=0.85;
fit = (a*double((beadsPerm)-1)).^1.5;%check
% Take mean over datablocks
meanR2Rosen = mean(r2Rosen, 2);
meanR2Perm = mean(r2Perm, 2);
% Plot loglog
fig = figure;
plot(beadsPerm, fit, '--k', ...
    beadsRosen, blockSizePerm(beadsRosen), 'og',...
    beadsRosen, meanR2Rosen(beadsRosen), '.r',...
    beadsPerm, meanR2Perm(beadsPerm), '.b');
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
ylabel('$R^2$','Interpreter','LaTex')%, 'FontSize',20);   % 'FontSize',20,'position',[-1.25,0.2]));
xlabel('N')
legend('fit','PERM', 'Standard rosenbluth')

%% 2nd plot
close all;
%beadsPerm=[min(beadsPerm,[],2):9, 10:5:40,50:10:90, 100:20:280, 300:50:max(beadsPerm,[],2)];
beadsPermNew = round(logspace(0,log10(double(max(beadsPerm,[],2))),50));  % Create the log-spaced index
beadsPermNew = unique(beadsPermNew);                   % Remove duplicate indices
sigmaPerm = std(r2Perm,1,2);
fit = (a*double((beadsPermNew)-1)).^1.5;
fig2 = figure;
plot(beadsPermNew, fit, '--k', ...
    beadsPermNew, blockSizePerm(beadsPermNew), 'og')
hold on;
errorbar(beadsPermNew,meanR2Perm(beadsPermNew),sigmaPerm(beadsPermNew))
ax = get(fig2,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
ylabel('$R^2$','Interpreter','LaTex')%, 'FontSize',20);   % 'FontSize',20,'position',[-1.25,0.2]));
xlabel('N')
legend('fit', 'Number of polymers','PERM')




