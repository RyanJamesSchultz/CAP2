% Simple script to demonstrate that k=1 order statistics are equivalent to the GR-MFD.
% See also Supplementary Text S2.
clear;

% Define some constants.
b=1.0;
m1=0.0;
m2=4;
k=1;
Nc=1;
d=1e-3;
YL=[1e-4 1e2];

% Derive some more variables.
m=m1:d:10;
a=log10(Nc);

% Expected distribution of M.
[Mpdf,Mcdf,Msvf]=GR_MFD(m, m1,m2, a,b, 'normalized'); Mpdf(m>m2)=0; Mcdf(m>m2)=1; Msvf(m>m2)=0;

% Expected distributions of Mlrg(i) & dMlrg(i).
[MLKpdf,MLKcdf,MLKsvf]=OSD(Mpdf,Mcdf,Msvf,k,Nc);
[dm,dMLKpdf,dMLKcdf,dMLKsvf]=dMD(MLKpdf,d);

% Simplified distribution of Mlrg(i) & dMlrg(i).
[MLKpdf2,MLKcdf2,MLKsvf2]=GR_MFD(m, m1,m2, a,Nc*b, 'normalized'); MLKpdf2(m>m2)=0; MLKcdf2(m>m2)=1; MLKsvf2(m>m2)=0;
[dMLKpdf2,dMLKcdf2,dMLKsvf2]=GR_MFD(dm, m1,m2, a,Nc*b, 'normalized'); dMLKpdf2(dm<0)=0; dMLKcdf2(dm<0)=1; dMLKsvf2(dm<0)=0;

% Plot.
figure(1); clf;

% Compare approaches for Mlrg.
subplot(2,1,1);
plot(m-m1,MLKpdf,'-k','DisplayName','Full PDF'); hold on;
plot(m-m1,MLKpdf2,'--r','DisplayName','Simple PDF');
plot(m-m1,100*MLKsvf,'-k','LineWidth',2,'DisplayName','Full SVF');
plot(m-m1,100*MLKsvf2,':r','LineWidth',2,'DisplayName','Simple SVF');
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL);
xlim([m1 inf])
legend();
title(['M_{LRG (1,',num2str(Nc),')}']);

% Compare approaches for dMlrg.
subplot(2,1,2);
plot(dm,dMLKpdf,'-k','DisplayName','Full PDF'); hold on;
plot(dm,dMLKpdf2,'--r','DisplayName','Simple PDF');
plot(dm,100*dMLKsvf,'-k','LineWidth',2,'DisplayName','Full SVF');
plot(dm,100*dMLKsvf2,':r','LineWidth',2,'DisplayName','Simple SVF');
xlabel('\DeltaM_{LRG}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
ylim(YL); 
xlim([0 inf])
legend();
title(['\DeltaM_{LRG (1,',num2str(Nc),')}']);

%
%figure(2); clf;
%semilogy(dm,dMLKpdf2./dMLKpdf)






