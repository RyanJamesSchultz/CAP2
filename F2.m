% Script to show composite metric distributions (i.e., M, dM, Mlrg, & dMlrg).
% Used to make Figure 2.
clear;

% Define some constants.
b=1.0;
m1=0.0;
m2=Inf;
Nc=100;
iL=[1 3 6 10 30 60 100];
d=1e-4;
YL=[1e-4 1e3];

% Derive some more variables.
kd=Nc-iL+1;
m=m1:d:10;
a=1e4;

% Expected distribution of M & dM.
[Mpdf,Mcdf,Msvf]=GR_MFD(m, m1,m2, a,b, 'normalized'); Mpdf(m>m2)=0; Mcdf(m>m2)=1; Msvf(m>m2)=0;
[dm,dMpdf,dMcdf,dMsvf]=dMD(Mpdf,d);

% Preallocate strucutres.
Sk=struct('MLpdf',[],'dMLpdf',[]);

% Loop over the k-index.
for k=1:length(kd)

    % Expected distributions of Mlrg(i,N) & dMlrg(i,N).
    [MLpdf,MLcdf,MLsvf]=OSD(Mpdf,Mcdf,Msvf,kd(k),Nc);
    [~,dMLpdf,dMLcdf,dMLsvf]=dMD(MLpdf,d);

    % Save the results.
    Sk(k).MLpdf=MLpdf;
    Sk(k).dMLpdf=dMLpdf;
end





% Plot.
figure(2); clf;

% Mlrg(i,N).
subplot(211);
plot(m-m1,Mpdf,'-k','LineWidth',2,'DisplayName','GR-MFD'); hold on;
for k=1:length(Sk)
    plot(m-m1,Sk(k).MLpdf,'DisplayName',['M_{LRG(',num2str(iL(k)),',',num2str(Nc),')}']);
end
xlabel('M-m_1'); ylabel('PDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('M_{LRG(i,N)}');
legend();

% dMlrg(i,N).
subplot(212);
plot(dm,dMpdf,'-k','LineWidth',2,'DisplayName','GR-MFD'); hold on;
for k=1:length(Sk)
    plot(dm,Sk(k).dMLpdf,'DisplayName',['\DeltaM_{LRG(',num2str(iL(k)),',',num2str(Nc),')}']);
end
xlabel('\DeltaM'); ylabel('PDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('\DeltaM_{LRG(i,N)}');
legend();
