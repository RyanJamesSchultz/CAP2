% Script to examine PDFs/SVFs of composite catalogue metrics.
% Used to make Figure S2.
clear;

% Define some constants.
b=1.0;
m1=0.0;
m2=Inf;
Nc=1e0;
d=1e-4;
kd=10;
YL=[1e-4 1e2];
ff=1;

% Derive some more variables.
m=m1:d:10;
a=log10(Nc);

% Expected distribution of M.
[Mpdf,Mcdf,Msvf]=GR_MFD(m, m1,m2, a,b, 'normalized'); Mpdf(m>m2)=0; Mcdf(m>m2)=1; Msvf(m>m2)=0;

% Expected distributions of Mlrg(i).
[ML1pdf,ML1cdf,ML1svf]=OSD(Mpdf,Mcdf,Msvf,1+ff,Nc+0+ff);
[ML2pdf,ML2cdf,ML2svf]=OSD(Mpdf,Mcdf,Msvf,1+ff,Nc+1+ff);
[ML3pdf,ML3cdf,ML3svf]=OSD(Mpdf,Mcdf,Msvf,1+ff,Nc+2+ff);
[ML4pdf,ML4cdf,ML4svf]=OSD(Mpdf,Mcdf,Msvf,1+ff,Nc+3+ff);
[MLHpdf,MLHcdf,MLHsvf]=OSD(Mpdf,Mcdf,Msvf,1+ff,Nc+(kd-1)+ff);

% Expected distributions of dM & dMlrg(i).
[dm,dMpdf,dMcdf,dMsvf]=dMD(Mpdf,d);
[~,dML1pdf,dML1cdf,dML1svf]=dMD(ML1pdf,d);
[~,dML2pdf,dML2cdf,dML2svf]=dMD(ML2pdf,d);
[~,dML3pdf,dML3cdf,dML3svf]=dMD(ML3pdf,d);
[~,dML4pdf,dML4cdf,dML4svf]=dMD(ML4pdf,d);
[~,dMLHpdf,dMLHcdf,dMLHsvf]=dMD(MLHpdf,d);




% Plot.
figure(52); clf;

% M
subplot(2,6,1);
plot(m-m1,Mpdf,'-k'); hold on;
plot(m-m1,100*Msvf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('M');
% Mlrg
subplot(2,6,2);
plot(m-m1,ML1pdf,'-k'); hold on;
plot(m-m1,100*ML1svf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('PDF/1-CDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['M_{LRG(',num2str(1+ff),',',num2str(Nc+0+ff),')}']);
% Mlrg(2)
subplot(2,6,3);
plot(m-m1,ML2pdf,'-k'); hold on;
plot(m-m1,100*ML2svf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('PDF/1-CDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['M_{LRG(',num2str(1+ff),',',num2str(Nc+1+ff),')}']);
% Mlrg(3)
subplot(2,6,4);
plot(m-m1,ML3pdf,'-k'); hold on;
plot(m-m1,100*ML3svf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('PDF/1-CDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['M_{LRG(',num2str(1+ff),',',num2str(Nc+2+ff),')}']);
% Mlrg(4)
subplot(2,6,5);
plot(m-m1,ML4pdf,'-k'); hold on;
plot(m-m1,100*ML4svf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('PDF/1-CDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['M_{LRG(',num2str(1+ff),',',num2str(Nc+3+ff),')}']);
% Mlrg(H)
subplot(2,6,6);
plot(m-m1,MLHpdf,'-k'); hold on;
plot(m-m1,100*MLHsvf,'-k','LineWidth',2);
xlabel('M-m_1'); ylabel('PDF/1-CDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['M_{LRG(',num2str(1+ff),',',num2str(Nc+kd+ff),')}']);

% dM
subplot(2,6,7);
plot(dm,dMpdf,'-k'); hold on;
plot(dm,100*dMsvf,'-k','LineWidth',2);
xlabel('\DeltaM'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('\DeltaM');
% dMlrg(1)
subplot(2,6,8);
plot(dm,dML1pdf,'-k'); hold on;
plot(dm,100*dML1svf,'-k','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+0+ff),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+0+ff),')}']);
% dMlrg(2)
subplot(2,6,9);
plot(dm,dML2pdf,'-k'); hold on;
plot(dm,100*dML2svf,'-k','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+1+ff),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+1+ff),')}']);
% dMlrg(3)
subplot(2,6,10);
plot(dm,dML3pdf,'-k'); hold on;
plot(dm,100*dML3svf,'-k','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+2+ff),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+2+ff),')}']);
% dMlrg(4)
subplot(2,6,11);
plot(dm,dML4pdf,'-k'); hold on;
plot(dm,100*dML4svf,'-k','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+3+ff),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+3+ff),')}']);
% dMlrg(H)
subplot(2,6,12);
plot(dm,dMLHpdf,'-k'); hold on;
plot(dm,100*dMLHsvf,'-k','LineWidth',2);
xlabel(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+kd+ff),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title(['\DeltaM_{LRG(',num2str(1+ff),',',num2str(Nc+kd+ff),')}']);

% Numerical error stuff.
figure(2); clf;
subplot(1,5,1);
semilogy(m-m1,abs(Mpdf-dMpdf(dm>=0))./Mpdf);
subplot(1,5,2);
semilogy(m-m1,abs(ML1pdf-dML1pdf(dm>=0))./ML1pdf);
subplot(1,5,3);
semilogy(m-m1,abs(ML2pdf-dML2pdf(dm>=0))./ML2pdf);
subplot(1,5,4);
semilogy(m-m1,abs(ML3pdf-dML3pdf(dm>=0))./ML3pdf);
subplot(1,5,5);
semilogy(m-m1,abs(ML4pdf-dML4pdf(dm>=0))./ML4pdf);





