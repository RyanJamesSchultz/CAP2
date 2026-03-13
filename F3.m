% Script to show composite metric distributions (i.e., M, dM, Mlrg, & dMlrg).
% Used to make Figure 3.
clear;

% Define some constants.
b=1.0;
m1=0.0;
m2=Inf;
Nc=[1 2 3 6 10 20 30 60 100];
k=1;
d=1e-4;
YL=[1e-4 1e3];

% Derive some more variables.
m=m1:d:10;
a=1e4;

% Expected distribution of M & dM.
[Mpdf,Mcdf,Msvf]=GR_MFD(m, m1,m2, a,b, 'normalized'); Mpdf(m>m2)=0; Mcdf(m>m2)=1; Msvf(m>m2)=0;
[dm,dMpdf,dMcdf,dMsvf]=dMD(Mpdf,d);

% Preallocate strucutres.
Sn=struct('MLpdf',[],'dMLpdf',[]);

% Loop over the k-index.
for n=1:length(Nc)
    
    % Expected distributions of Mlrg(i,N) & dMlrg(i,N).
    [MLpdf,MLcdf,MLsvf]=OSD(Mpdf,Mcdf,Msvf,k,Nc(n));
    [~,dMLpdf,dMLcdf,dMLsvf]=dMD(MLpdf,d);
    
    % Save the results.
    Sn(n).MLpdf=MLpdf;
    Sn(n).dMLpdf=dMLpdf;
end





% Plot.
figure(3); clf;

% Mlrg(i,N).
subplot(211);
plot(m-m1,Mpdf,':k','LineWidth',2,'DisplayName','GR-MFD'); hold on;
for n=1:length(Sn)
    plot(m-m1,Sn(n).MLpdf,'DisplayName',['M_{LRG(',num2str(Nc(n)-k+1),',',num2str(Nc(n)),')}']);
end
xlabel('M-m_1'); ylabel('PDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('M_{LRG(i,N)}');
legend();

% dMlrg(i,N).
subplot(212);
plot(dm,dMpdf,':k','LineWidth',2,'DisplayName','GR-MFD'); hold on;
for n=1:length(Sn)
    plot(dm,Sn(n).dMLpdf,'DisplayName',['\DeltaM_{LRG(',num2str(Nc(n)-k+1),',',num2str(Nc(n)),')}']);
end
xlabel('\DeltaM'); ylabel('PDF');
set(gca,'YScale','log');
ylim(YL);
xlim([0 4]);
title('\DeltaM_{LRG(i,N)}');
legend();
