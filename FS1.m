% Script to compare a numerical distribution of dM versus the analytical expectation.
% Used to make Figure S1.
clear;

% Define some constants.
Nc=1e6;
m1=0.0;
b=1.0;
mr=-(0.5:0.1:3.0);

% Derive some more variables.
m2=linspace(3,3,Nc);
a=log10(Nc);
m=m1:0.02:10;
dm=0:0.02:10;

% Get a catalogue.
M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);

% Get the expected distributions.
[pdf,cdf,svf]=GR_MFD(m,m1,max(m2),a,b,'norm');
[PDF,CDF,SVF]=GR_dM0(dm,max(m2)-m1,a,b,'norm');

% Get the two sets of magnitude jumps.
dM=diff([m1,M]);
dMp=dM(dM>0);
dMn=abs(dM(dM<0));
dMa=abs(dM);




% Plot.
figure(51); clf;
% The GR-MFD.
subplot(311);
histogram(M-m1,m-m1,'Normalization','pdf'); hold on;
plot(m-m1,pdf,':k');
plot(m-m1,100*svf,':k','LineWidth',2);
plot(sort(M-m1),100*getESVF(M),'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 max(m2-m1)+1]); YL=ylim();
title('GR-MFD');
% Positive jumps.
subplot(312);
histogram(dMp,dm,'Normalization','pdf'); hold on;
plot(m-m1,pdf,':k');
plot(dm,PDF,'-k');
plot(m-m1,100*svf,':k','LineWidth',2);
plot(dm,100*SVF,'-k','LineWidth',2);
plot(sort(dMp),100*getESVF(dMp),'-b','LineWidth',2);
xlabel('\DeltaM'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 max(m2-m1)+1]); YL=ylim();
title('\DeltaM(+)');
% Negative jumps.
subplot(313);
histogram(dMn,dm,'Normalization','pdf'); hold on;
plot(m-m1,pdf,':k');
plot(dm,PDF,'-k');
plot(m-m1,100*svf,':k','LineWidth',2);
plot(dm,100*SVF,'-k','LineWidth',2);
plot(sort(dMn),100*getESVF(dMn),'-b','LineWidth',2);
xlabel('\DeltaM'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 max(m2-m1)+1]); YL=ylim();
title('\DeltaM(-)');





%%%% SUBROUNTINES.

% Given a sample dataset, find the empirical SVF.
function [eSVF]=getESVF(x)
  eSVF=(1-(1:length(x))/length(x));
end