% Script to find the order statistic summation factors for the empirical 
% distribution of dMlrg(i).
% This is important for the MLE fitting of Wj in Section 2.5.3.
clear;
rng('shuffle');

% Define some constants.
Nc=1e4;
Ns=1e4;
m1=0.0;
m2=Inf;
b=1.0;
d1=1e-4;
d2=5e-3;
kd=20;

% Derive some more variables.
a=log10(Nc);
T=1:Nc;

% Preallocate.
dM_lrgK=[];

% Loop.
disp('--- Loop');
for i=1:Ns
    
    % Draw a catalogue.
    Mb=GR_MFD_Rand(m1,m2 , a,b, [1 Nc]);
    
    % Get the biggest event sequence (for the bounded catalogue).
    MlrgK_i=OrderStatistic(Mb,length(Mb)-(kd-1),'unique'); MlrgK_i=[m1,MlrgK_i]; MlrgK_i(isnan(MlrgK_i))=[];
    dM_lrgK=[dM_lrgK, diff(MlrgK_i)];
    
end

% Make plotting axes. 
mrB1=min(Mb):d1:min([m2,10]);
mrB2=min(Mb):d2:min([m2,10]);

% Fit weights for best dMlrg PDF.
disp('--- Fit');
wg=1./(1:(kd+1)).^(1.1); wg=wg/wg(end); wg(end-5:end)=[];

% Fit norm type.
type='norm';
nLL = @(w1)  -GR_NLE_LL(mrB1,dM_lrgK, m1,max(m2), a,b, kd,type,[w1,0,0,0,0,0,1]);
%w = fmincon(nLL,wg,[],[],[],[],1e-4*ones([1 k-5]),1e4*ones([1 k-5]));
%w=[w,0,0,0,0,0,1];
% Fit fudge type.
%type='fudge';
%nLL = @(w)  -GR_NLE_LL(mrB1,dM_lrgK, m1,max(m2), a,b, k,type,w);
%w = fmincon(nLL,1,[],[],[],[],1e-1,1e4);

% k=2
%w=[2.2 0.9 1];
%w=[1.55 0 1];
%w=[1 0.00 0.70];

% k=3
%w=[8 0 7 1];
%w=[4.6 0.9 3.2 1];
%w=[1.75 1.5 0 1];
%w=[1 0 1 0];

% k=4
%w=[1.61 2.01 0 0 1];
%w=[1.5 1.5 0 0 1];
%w=[1.19    0.73    0    1         0];
%w=[4.98    0    5.04    0    1];
%w=[2.13    1.36    0.84    0    1];
%w=[4.70 4.37 0 3.78 0];
%w=[3.92    0    3.69    1    0];

% k=5
%w=[4.94    0    5.66    0    0    1];
%w=[1.5 2.6 0 0 0 1];
%w=[5.43    0.33    6.29    0    0    1];
%w=[1.93    1.96    0    0.65    1    0];
%w=[1.29 1.83 0 0 1 0];

% k=6
%w=[3.40 1.11 3.47 0 0    0    1];
%w=[2.98 1.34 2.44 0 0.01 1.22 0];
%w=[0.95    0.08    1.12    0    0    0.24         0];
%w=[488   42  578    0    0  123    1];
%w=[547 660 0.4 116 138 218 1];
%w=[26.44 29.66 5.40 1.66 5.04 12.45 0];
%w=[12 6.65 9.11 0.06 2.64 3.04 0];
%w=[5.14    0    6.95    0    0    0    1];
%w=[4.12    0    5.30    0    0    1    0];

% k=7
%w=[5.50    0    7.83    0    0    0    0    1];
%w=[5.27    0    7.72    0    0    0    0    1];
%w=[4.37 0.69    5.78    0    0    0    0    1];
%w=[4.60 4.87    0    3.92    0    0    0    1];

% k=8
%w=[5.79 6.57    0    5.90    0    0    0    0    1];
%w=[4.43 5.34    0    4.12    0    0    0    0    1];
%w=[4.48 0.80    6.62    0    0    0    0    0    1];
%w=[5.35    0    8.11    0    0    0    0    0    1];

% k=9
%w=[10.93 8.34 0.00 13.08 0 0 0 0 0 1];
%w=[ 6.75 0.00 9.18  1.61 0 0 0 0 0 1];
%w=[ 5.91 5.85 0.00  6.41 0 0 0 0 0 1]; %biased.
%w=[ 5.68 0.00 7.35  1.15 0 0 0 0 0 1]; %biased.
%w=[ 5.98 0.00 7.33  1.59 0 0 0 0 0 1]; %biased.

% k=10
w=[48.31 65.80  0.53 0.70 54.45 0 0 0 0 0 1];
w=[10.14  0.00 10.38 5.79  0.00 0 0 0 0 0 1];
w=[10.32  1.71 12.17 0.00  4.20 0 0 0 0 0 1];
w=[ 8.12  0.00  9.20 4.07  0.00 0 0 0 0 0 1];
%w=[11.85  0.00 10.59 6.72  1.03 0 0 0 0 0 1]; %biased.


% Expected M & dMlrg(i) PDFs.
[MBpdf,MBcdf,MBsvf]=GR_MFD(mrB1, m1,max(m2), a,b,'normalized');
[dmB,dMLKBpdf,dMLKBcdf,dMLKBsvf]=GR_NLE(mrB1, m1,max(m2), a,b, kd,type,w);

% Empirical distribution functions.
[dMecdf,dMesvf]=getESVF(dM_lrgK);
[Mecdf,Mesvf]=getESVF(Mb);

% Quick KS-test.
rng(3);
[~,I]=unique(dMLKBcdf);
dMLs=interp1(dMLKBsvf(I),dmB(I),rand([1 Nc]),'linear');
[~,KSp]=kstest2(dMLs,dM_lrgK);

%Report values.
w
KSp



% Plot.
figure(1); clf;
Mp=max(dM_lrgK);

% GR-MFD.
subplot(2,1,1);
histogram(Mb-m1,mrB2-m1,'Normalization','pdf'); hold on;
plot(mrB1-m1,MBpdf,'-k');
plot(mrB1-m1,100*MBsvf,'-k','LineWidth',2);
plot(sort(Mb-m1),100*Mesvf,'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 Mp]); ylim([interp1(mrB1-m1,MBpdf,Mp,'linear'),1.3*max(100*Mesvf)]);
title('M');
% dMlrg,k.
subplot(2,1,2);
histogram(dM_lrgK,mrB2-m1,'Normalization','pdf'); hold on;
plot(dmB,dMLKBpdf,'-k');
plot(dmB,100*dMLKBsvf,'-k','LineWidth',2);
plot(sort(dM_lrgK),100*dMesvf,'-b','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(kd),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 Mp]); ylim([interp1(dmB,dMLKBpdf,Mp,'linear'),1.3*max(100*dMesvf)]);
title(['\DeltaM_{LRG(',num2str(kd),')}']);

% Residual.
figure(2); clf;
subplot(2,1,1);
plot([0 max(dM_lrgK)+0.1],[0 0],'-k','LineWidth',2); hold on;
plot(sort(dM_lrgK),100*(dMesvf-interp1(dmB,    dMLKBsvf,sort(dM_lrgK),'linear'))./dMesvf,'-b','LineWidth',2);
plot(sort(Mb-m1),  100*( Mesvf-interp1(mrB1-m1,MBsvf,   sort(Mb-m1),   'linear'))./Mesvf,':b');
%plot(sort(dM_lrgK),100*(dMecdf-interp1(dmB,    dMLKBcdf,sort(dM_lrgK),'linear'))./dMecdf,'-r','LineWidth',2);
%plot(sort(Mb-m1),  100*( Mecdf-interp1(mrB1-m1,MBcdf,   sort(Mb-m1),   'linear'))./Mecdf,':r');
xlabel(['\DeltaM_{LRG(',num2str(kd),')}']); ylabel('Residual Error (%)');
ylim([-15 +15]); 
xlim([0 max(dM_lrgK)+0.1])
subplot(2,1,2);
plot([0 max(dM_lrgK)+0.1],[0 0],'-k'); hold on;
plot(sort(dM_lrgK),abs(100*(dMesvf-interp1(dmB,    dMLKBsvf,sort(dM_lrgK),'linear'))./dMesvf),'-b','LineWidth',2);
plot(sort(Mb-m1),  abs(100*( Mesvf-interp1(mrB1-m1,MBsvf,   sort(Mb-m1),   'linear'))./Mesvf),':b');
%plot(sort(dM_lrgK),100*(dMecdf-interp1(dmB,    dMLKBcdf,sort(dM_lrgK),'linear'))./dMecdf,'-r','LineWidth',2);
%plot(sort(Mb-m1),  100*( Mecdf-interp1(mrB1-m1,MBcdf,   sort(Mb-m1),   'linear'))./Mecdf,':r');
xlabel(['\DeltaM_{LRG(',num2str(kd),')}']); ylabel('Residual Error (%)');
set(gca,'YScale','log');
xlim([0 max(dM_lrgK)+0.1])





%%%% SUBROUNTINES.

% Given a sample dataset, find the empirical SVF.
function [eCDF, eSVF]=getESVF(x)
  eCDF=(1:length(x))/length(x);
  eSVF=1-eCDF;
end

