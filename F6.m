% Plots the test bench MLE-fitting results.
% See also script_MakeMLE_TestBench.m.
% Used to make Figures 6 & S8.
clear;

% Load in the data from the test bench.
load('TestBenchMLE_B10.mat','S'); % Composite form of nLL.
%load('TestBenchMLE_C10.mat','S'); % Best/combined form of nLL.
%load('TestBenchMLE_A10.mat', 'S'); % Approximate nLL.

% Get results.
Mmax_Bavg1=arrayfun(@(S) prctile(abs(S.Mmax_bis),50),S);
Mmax_Bavg2=arrayfun(@(S) abs(prctile(S.Mmax_bis,50)),S);
Mmax_Bp10=arrayfun(@(S) prctile(abs(S.Mmax_bis),10),S);
Mmax_Bp90=arrayfun(@(S) prctile(abs(S.Mmax_bis),90),S);
Mmax_Bstd=arrayfun(@(S) std(S.Mmax_bis),S);
Mmax_Savg=arrayfun(@(S) prctile(S.Mmax_std,50),S);
Mmax_Sp10=arrayfun(@(S) prctile(S.Mmax_std,10),S);
Mmax_Sp90=arrayfun(@(S) prctile(S.Mmax_std,90),S);
dMexp=arrayfun(@(S) S.dMexp,S); dMexp(isinf(dMexp))=1;
Nc=arrayfun(@(S) S.Nc,S);
Ns=length(S(1,1).Mmax_bis);

% Report test bench values.
Mmax_Bavg1
Mmax_Bavg2
Mmax_Bstd
Mmax_Savg
%dMexp
%Nc




% Plot.
colours={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#A2142F'};

% Plot variance/bias trends in MLE-fitted Mmax.
figure(6); clf;
% Bias in Mmax estimates.
subplot(211);
semilogy([min(dMexp(1,:)) max(dMexp(1,:))],[0 0],'-k','HandleVisibility','off'); hold on;
for i=1:size(Nc,1)
    semilogy(dMexp(i,:),Mmax_Bavg1(i,:),'-','Color',colours{i},'DisplayName',['MAD Estimate (10^{',num2str(log10(Nc(i,1))),'})']);
    semilogy(dMexp(i,:),Mmax_Bp10(i,:),':','Color',colours{i},'HandleVisibility','off');
    semilogy(dMexp(i,:),Mmax_Bp90(i,:),':','Color',colours{i},'HandleVisibility','off');
    semilogy(dMexp(i,:),Mmax_Bavg2(i,:),'--','Color',colours{i},'DisplayName',['AMD Estimate (10^{',num2str(log10(Nc(i,1))),'})']);
end
xlabel('Degree of Truncation, dMx'); ylabel('M_{MAX} Absolute Bias');
xlim([-1 +1]); ylim([5e-3 2e0]);
legend('Location','northwest');
% Variance in Mmax estimates.
subplot(212);
semilogy([min(dMexp(1,:)) max(dMexp(1,:))],[0 0],'-k','HandleVisibility','off'); hold on;
for i=1:size(Nc,1)
    semilogy(dMexp(i,:),Mmax_Savg(i,:),'-','Color',colours{i},'DisplayName',['Estimate (10^{',num2str(log10(Nc(i,1))),'})']);
    semilogy(dMexp(i,:),Mmax_Sp10(i,:),':','Color',colours{i},'HandleVisibility','off');
    semilogy(dMexp(i,:),Mmax_Sp90(i,:),':','Color',colours{i},'HandleVisibility','off');
    semilogy(dMexp(i,:),Mmax_Bstd(i,:),'--','Color',colours{i},'DisplayName',['True (10^{',num2str(log10(Nc(i,1))),'})']);
end
xlabel('Degree-of-Truncation, \deltaM'); ylabel('M_{MAX} Standard Deviation');
xlim([-1 +1]); ylim([5e-3 2e0]);
legend('Location','northwest');

% Bonus plots of individual MLE-fit variance/bias metrics.
figure(1); clf;
i=1; j=3;
% Historgrams of Mmax absolute bias.
subplot(511);
histogram(abs(S(i,j).Mmax_bis),round(sqrt(5*Ns))); hold on;
plot([0 0],ylim(),'-k');
plot(Mmax_Bavg1(i,j)*[1 1],ylim(),'--b');
xlabel('M_{MAX} Bias'); ylabel('Counts');
% Historgrams of Mmax bias.
subplot(512);
histogram((S(i,j).Mmax_bis),round(sqrt(5*Ns))); hold on;
plot([0 0],ylim(),'-k');
plot(median(S(i,j).Mmax_bis)*[1 1],ylim(),'--b');
xlabel('M_{MAX} Bias'); ylabel('Counts');
% Histograms of Mmax standard deviation.
subplot(513);
histogram(S(i,j).Mmax_std,round(sqrt(5*Ns))); hold on;
plot(Mmax_Bstd(i,j)*[1 1],ylim(),'-k');
plot(Mmax_Savg(i,j)*[1 1],ylim(),'--b');
xlabel('M_{MAX} Standard Deviation'); ylabel('Counts');
subplot(5,1,[4 5]);
mx=max([max(abs(S(i,j).Mmax_bis)),max(S(i,j).Mmax_std)]);
plot([0 mx],[0 mx],'-k'); hold on;
plot(abs(S(i,j).Mmax_bis),S(i,j).Mmax_std,'ob');
xlabel('M_{MAX} Absolute Bias'); ylabel('M_{MAX} Standard Deviation');
% Report plotted values.
Nc(i,j)
dMexp(i,j)
