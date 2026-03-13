% Script to examine the dMlrg(i) distributions numerically and synthetically.
% Used to make Figure S3.
clear;

% Define some constants.
Nc=1e4;
Ns=1e3;
m1=0.0;
b=1.0;
d1=1e-4;
d2=0.1;
kd=10;
type='norm';
O=-1;

% Derive some more variables.
m2=linspace(2.5,4.0,Nc);
a=log10(Nc);
T=1:Nc;

% Preallocate.
Mlrg1U=[]; Mlrg2U=[]; Mlrg3U=[]; Mlrg4U=[]; MlrgHU=[];
Mlrg1B=[]; Mlrg2B=[]; Mlrg3B=[]; Mlrg4B=[]; MlrgHB=[];

dM_lrg1U=[]; dM_lrg2U=[]; dM_lrg3U=[]; dM_lrg4U=[]; dM_lrgHU=[];
dM_lrg1B=[]; dM_lrg2B=[]; dM_lrg3B=[]; dM_lrg4B=[]; dM_lrgHB=[];

% Loop.
for i=1:Ns
    
    % Draw a catalogue.
    Mb=GR_MFD_Rand(m1,m2 , a,b, [1 Nc]);
    Mu=GR_MFD_Rand(m1,Inf, a,b, [1 Nc]);
    
    % Get the dMlrg(i) seqeunce (for the bounded catalogue).
    Mlrg1_i=OrderStatistic(Mb,length(Mb)-0, 'unique'); Mlrg1_i=[m1,Mlrg1_i]; Mlrg1_i(isnan(Mlrg1_i))=[];
    Mlrg2_i=OrderStatistic(Mb,length(Mb)-1, 'unique'); Mlrg2_i=[m1,Mlrg2_i]; Mlrg2_i(isnan(Mlrg2_i))=[];
    Mlrg3_i=OrderStatistic(Mb,length(Mb)-2, 'unique'); Mlrg3_i=[m1,Mlrg3_i]; Mlrg3_i(isnan(Mlrg3_i))=[];
    Mlrg4_i=OrderStatistic(Mb,length(Mb)-3, 'unique'); Mlrg4_i=[m1,Mlrg4_i]; Mlrg4_i(isnan(Mlrg4_i))=[];
    MlrgH_i=OrderStatistic(Mb,length(Mb)-(kd-1),'unique'); MlrgH_i=[m1,MlrgH_i]; MlrgH_i(isnan(MlrgH_i))=[];
    Mlrg1B=[Mlrg1B, Mlrg1_i];
    Mlrg2B=[Mlrg2B, Mlrg2_i];
    Mlrg3B=[Mlrg3B, Mlrg3_i];
    Mlrg4B=[Mlrg4B, Mlrg4_i];
    MlrgHB=[MlrgHB, MlrgH_i];
    dM_lrg1B=[dM_lrg1B, diff(Mlrg1_i)];
    dM_lrg2B=[dM_lrg2B, diff(Mlrg2_i)];
    dM_lrg3B=[dM_lrg3B, diff(Mlrg3_i)];
    dM_lrg4B=[dM_lrg4B, diff(Mlrg4_i)];
    dM_lrgHB=[dM_lrgHB, diff(MlrgH_i)];
    
    % Get the dMlrg(i) seqeunce (for the unbounded catalogue).
    Mlrg1_i=OrderStatistic(Mu,length(Mu)-0, 'unique'); Mlrg1_i=[m1,Mlrg1_i]; Mlrg1_i(isnan(Mlrg1_i))=[];
    Mlrg2_i=OrderStatistic(Mu,length(Mu)-1, 'unique'); Mlrg2_i=[m1,Mlrg2_i]; Mlrg2_i(isnan(Mlrg2_i))=[];
    Mlrg3_i=OrderStatistic(Mu,length(Mu)-2, 'unique'); Mlrg3_i=[m1,Mlrg3_i]; Mlrg3_i(isnan(Mlrg3_i))=[];
    Mlrg4_i=OrderStatistic(Mu,length(Mu)-3, 'unique'); Mlrg4_i=[m1,Mlrg4_i]; Mlrg4_i(isnan(Mlrg4_i))=[];
    MlrgH_i=OrderStatistic(Mu,length(Mu)-(kd-1),'unique'); MlrgH_i=[m1,MlrgH_i]; MlrgH_i(isnan(MlrgH_i))=[];
    Mlrg1U=[Mlrg1U, Mlrg1_i];
    Mlrg2U=[Mlrg2U, Mlrg2_i];
    Mlrg3U=[Mlrg3U, Mlrg3_i];
    Mlrg4U=[Mlrg4U, Mlrg4_i];
    MlrgHU=[MlrgHU, MlrgH_i];
    dM_lrg1U=[dM_lrg1U, diff(Mlrg1_i)];
    dM_lrg2U=[dM_lrg2U, diff(Mlrg2_i)];
    dM_lrg3U=[dM_lrg3U, diff(Mlrg3_i)];
    dM_lrg4U=[dM_lrg4U, diff(Mlrg4_i)];
    dM_lrgHU=[dM_lrgHU, diff(MlrgH_i)];
    
end

% Test the sample distributions for similarity.
[~,Pu ]=kstest2(dM_lrg1U ,Mu-m1) % Indistinguishable.
[~,Pb ]=kstest2(dM_lrg1B ,Mb-m1) % Distinguishable.

% Expected M & Mlrg PDFs.
mrU1=min(Mb):d1:10;
mrB1=min(Mb):d1:max(m2);
mrU2=min(Mb):d2:max([Mu dM_lrg1U]);
mrB2=min(Mb):d2:max(m2);
[MUpdf,MUcdf,MUsvf]=GR_MFD(mrU1, m1,    Inf, a,b,'normalized');
[MBpdf,MBcdf,MBsvf]=GR_MFD(mrB1, m1,max(m2), a,b,'normalized');
[ML1Updf,ML1Ucdf,ML1Usvf]=OSD(MUpdf,MUcdf,MUsvf,1,1);
[ML1Bpdf,ML1Bcdf,ML1Bsvf]=OSD(MBpdf,MBcdf,MBsvf,1,1);
[ML2Updf,ML2Ucdf,ML2Usvf]=OSD(MUpdf,MUcdf,MUsvf,1.5,2.5);
[ML2Bpdf,ML2Bcdf,ML2Bsvf]=OSD(MBpdf,MBcdf,MBsvf,1.5,2.5);
[ML3Updf,ML3Ucdf,ML3Usvf]=OSD(MUpdf,MUcdf,MUsvf,2,4);
[ML3Bpdf,ML3Bcdf,ML3Bsvf]=OSD(MBpdf,MBcdf,MBsvf,2,4);
[ML4Updf,ML4Ucdf,ML4Usvf]=OSD(MUpdf,MUcdf,MUsvf,2,5);
[ML4Bpdf,ML4Bcdf,ML4Bsvf]=OSD(MBpdf,MBcdf,MBsvf,2,5);
[MLHUpdf,MLHUcdf,MLHUsvf]=OSD(MUpdf,MUcdf,MUsvf,1+1,kd+1);
[MLHBpdf,MLHBcdf,MLHBsvf]=OSD(MBpdf,MBcdf,MBsvf,1+1,kd+1);

% Expected distributions of dM & dMlrg(i).
[dmU,dMUpdf,dMUcdf,dMUsvf]=dMD(MUpdf,d1);
[dmB,dMBpdf,dMBcdf,dMBsvf]=dMD(MBpdf,d1);
[dMBL,dML1Updf,dML1Ucdf,dML1Usvf]=GR_NLE(mrU1, m1,    Inf, a,b,  1,type,O);
[   ~,dML1Bpdf,dML1Bcdf,dML1Bsvf]=GR_NLE(mrU1, m1,max(m2), a,b,  1,type,O);
[   ~,dML2Updf,dML2Ucdf,dML2Usvf]=GR_NLE(mrU1, m1,    Inf, a,b,  2,type,O);
[   ~,dML2Bpdf,dML2Bcdf,dML2Bsvf]=GR_NLE(mrU1, m1,max(m2), a,b,  2,type,O);
[   ~,dML3Updf,dML3Ucdf,dML3Usvf]=GR_NLE(mrU1, m1,    Inf, a,b,  3,type,O);
[   ~,dML3Bpdf,dML3Bcdf,dML3Bsvf]=GR_NLE(mrU1, m1,max(m2), a,b,  3,type,O);
[   ~,dML4Updf,dML4Ucdf,dML4Usvf]=GR_NLE(mrU1, m1,    Inf, a,b,  4,type,O);
[   ~,dML4Bpdf,dML4Bcdf,dML4Bsvf]=GR_NLE(mrU1, m1,max(m2), a,b,  4,type,O);
[   ~,dMLHUpdf,dMLHUcdf,dMLHUsvf]=GR_NLE(mrU1, m1,    Inf, a,b, kd,type,O);
[   ~,dMLHBpdf,dMLHBcdf,dMLHBsvf]=GR_NLE(mrU1, m1,max(m2), a,b, kd,type,O);




% Plot.
figure(53); clf;

% Regular unbounded GR-MFD.
subplot(2,6,1);
histogram(Mu-m1,mrU2-m1,'Normalization','pdf'); hold on;
plot(mrU1-m1,MUpdf,'-k');
plot(mrU1-m1,100*MUsvf,'-k','LineWidth',2);
plot(sort(Mu-m1),100*getESVF(Mu),'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); YL=ylim();
title('Unbounded M');
% Unbounded dMlrg.
subplot(2,6,2);
histogram(dM_lrg1U,mrU2-m1,'Normalization','pdf'); hold on;
plot(dmU,dML1Updf,'-k');
plot(dmU,100*dML1Usvf,'-k','LineWidth',2);
plot(sort(dM_lrg1U),100*getESVF(dM_lrg1U),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Unbounded \DeltaM_{LRG}');
% Unbounded dMlrg(2).
subplot(2,6,3);
histogram(dM_lrg2U,mrU2-m1,'Normalization','pdf'); hold on;
plot(dmU,dML2Updf,'-k');
plot(dmU,100*dML2Usvf,'-k','LineWidth',2);
plot(sort(dM_lrg2U),100*getESVF(dM_lrg2U),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(2)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Unbounded \DeltaM_{LRG(2)}');
% Unbounded dMlrg(3).
subplot(2,6,4);
histogram(dM_lrg3U,mrU2-m1,'Normalization','pdf'); hold on;
plot(dmU,dML3Updf,'-k');
plot(dmU,100*dML3Usvf,'-k','LineWidth',2);
plot(sort(dM_lrg3U),100*getESVF(dM_lrg3U),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(3)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Unbounded \DeltaM_{LRG(3)}');
% Unbounded dMlrg(4).
subplot(2,6,5);
histogram(dM_lrg4U,mrU2-m1,'Normalization','pdf'); hold on;
plot(dmU,dML4Updf,'-k');
plot(dmU,100*dML4Usvf,'-k','LineWidth',2);
plot(sort(dM_lrg4U),100*getESVF(dM_lrg4U),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(4)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Unbounded \DeltaM_{LRG(4)}');
% Unbounded dMlrg(H).
subplot(2,6,6);
histogram(dM_lrgHU,mrU2-m1,'Normalization','pdf'); hold on;
plot(dmU,dMLHUpdf,'-k');
plot(dmU,100*dMLHUsvf,'-k','LineWidth',2);
plot(sort(dM_lrgHU),100*getESVF(dM_lrgHU),'-b','LineWidth',2);
xlabel(['M_{LRG(',num2str(kd),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 1]); ylim(YL);
title(['Unbounded \DeltaM_{LRG(',num2str(kd),')}']);

% Bounded GR-MFD.
subplot(2,6,7);
histogram(Mb-m1,mrU2-m1,'Normalization','pdf'); hold on;
plot(mrB1-m1,MBpdf,'-k');
plot(mrB1-m1,100*MBsvf,'-k','LineWidth',2);
plot(sort(Mb-m1),100*getESVF(Mb),'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
plot((max(m2)-m1)*[1 1],ylim(),':k');
xlim([0 4]); ylim(YL);
title('Bounded M');
% Bounded dMlrg.
subplot(2,6,8);
histogram(dM_lrg1B,mrU2-m1,'Normalization','pdf'); hold on;
plot(dMBL,dML1Bpdf,'-k');
plot(dMBL,100*dML1Bsvf,'-k','LineWidth',2);
plot(sort(dM_lrg1B),100*getESVF(dM_lrg1B),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Bounded \DeltaM_{LRG}');
% Bounded dMlrg(2).
subplot(2,6,9);
histogram(dM_lrg2B,mrU2-m1,'Normalization','pdf'); hold on;
plot(dMBL,dML2Bpdf,'-k');
plot(dMBL,100*dML2Bsvf,'-k','LineWidth',2);
plot(sort(dM_lrg2B),100*getESVF(dM_lrg2B),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(2)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Bounded \DeltaM_{LRG(2)}');
% Bounded dMlrg(3).
subplot(2,6,10);
histogram(dM_lrg3B,mrU2-m1,'Normalization','pdf'); hold on;
plot(dMBL,dML3Bpdf,'-k');
plot(dMBL,100*dML3Bsvf,'-k','LineWidth',2);
plot(sort(dM_lrg3B),100*getESVF(dM_lrg3B),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(3)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Bounded \DeltaM_{LRG(3)}');
% Bounded dMlrg(4).
subplot(2,6,11);
histogram(dM_lrg4B,mrU2-m1,'Normalization','pdf'); hold on;
plot(dMBL,dML4Bpdf,'-k');
plot(dMBL,100*dML4Bsvf,'-k','LineWidth',2);
plot(sort(dM_lrg4B),100*getESVF(dM_lrg4B),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG(4)}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
title('Bounded \DeltaM_{LRG(4)}');
% Bounded dMlrg(H).
subplot(2,6,12);
histogram(dM_lrgHB,mrU2-m1,'Normalization','pdf'); hold on;
plot(dMBL,dMLHBpdf,'-k');
plot(dMBL,100*dMLHBsvf,'-k','LineWidth',2);
plot(sort(dM_lrgHB),100*getESVF(dM_lrgHB),'-b','LineWidth',2);
xlabel(['M_{LRG(',num2str(kd),')}']); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 1]); ylim(YL);
title(['Unbounded \DeltaM_{LRG(',num2str(kd),')}']);






%%%% SUBROUNTINES.

% Given a sample dataset, find the empirical SVF.
function [eSVF]=getESVF(x)
  eSVF=(1-(1:length(x))/length(x));
end
