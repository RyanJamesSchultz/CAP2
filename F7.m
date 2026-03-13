% Display an example of the EW-test.
% Used to make Figures 7, 8, & S9.
clear;
rng(2);

% Define some constants.
Nc=1e4;
Nr=0;
m1=0.0;
b=1.0;
kd=10;
reshuffle_flag='none';
smooth_flag='mean';
GRtype_flag='M';
NLEtype_flag='comp';

% Derive some more variables.
m2=linspace(2,4,Nc);
a=log10(Nc);
Mmax=max(m2);

% Load in lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');
%Sa=S;

% Get a catalogue.
M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);

% Define these order statistics, so I can see their values.
D=OrderStatisticsI(M,kd,'unique');
for j=1:length(D)
    D(j).Mlrg=[m1,D(j).Mlrg];
    D(j).Mlrg(isnan(D(j).Mlrg))=[];
    D(j).dMlrg=diff(D(j).Mlrg);
    D(j).Mlrg(end)=[];
end

% Prep for Mmax guesses.
Kgr=2+1;
i=1; name={};

% Make the Mmax guesses.
name{i}='MA1 (Linear)';
Sm(i).Mmax=m2;
Sm(i).K=Kgr+2;
i=i+1;
name{i}='MA2 (Linear + error)';
Sm(i).Mmax=m2+0.1;
Sm(i).K=Kgr+2;
i=i+1;
name{i}='MB1 (Log-Linear)';
c0=0.7;
c1=(max(m2)-c0)/log(Nc);
Sm(i).Mmax=c1*log(1:Nc)+c0;
Sm(i).K=Kgr+2;
i=i+1;
name{i}='MB2 (Log-Linear + error)';
Sm(i).Mmax=Sm(i-1).Mmax+0.1;
Sm(i).K=Kgr+2;
i=i+1;
name{i}='MC1 (Constant M_{LRG})';
Sm(i).Mmax=max(M)*ones(size(M))+0.1;
Sm(i).K=Kgr+1;
i=i+1;
name{i}='MC2 (Constant M_{LRG} + error)';
Sm(i).Mmax=max(M)+0.2;
Sm(i).K=Kgr+1;
i=i+1;
name{i}='MD1 (M_{LRG} Envelope)';
Sm(i).Mmax=cummax(M)+0.1;
Sm(i).K=Kgr+2*length(D(1).I);
i=i+1;
name{i}='MD2 (M_{LRG} Envelope + error)';
Sm(i).Mmax=cummax(M)+0.2;
Sm(i).K=Kgr+2*length(D(1).I);
i=i+1;
name{i}='Unbound';
Sm(i).Mmax=Inf;
Sm(i).K=Kgr+0;
i=i+1;

% Do the EW-tests.
tic;
W10=EnsembleW(M,m1,Sm,b,kd,Nr,reshuffle_flag,smooth_flag, GRtype_flag,NLEtype_flag,Sa,ma,ba,m2a);
toc;
W0=EnsembleW(M,m1,Sm,b,kd,Nr,reshuffle_flag,smooth_flag, GRtype_flag,'none',      Sa,ma,ba,m2a);
W1=EnsembleW(M,m1,Sm,b,1 ,Nr,reshuffle_flag,smooth_flag, GRtype_flag,NLEtype_flag,Sa,ma,ba,m2a);
W5=EnsembleW(M,m1,Sm,b,5 ,Nr,reshuffle_flag,smooth_flag, GRtype_flag,NLEtype_flag,Sa,ma,ba,m2a);

% Make weight matrices.
Wb10=[]; Wb5=[]; Wb1=[]; Wb0=[];
for j=1:length(W10)
    Wb10=[Wb10;W10(j).W];
    Wb5 =[ Wb5; W5(j).W];
    Wb1 =[ Wb1; W1(j).W];
    Wb0 =[ Wb0; W0(j).W];
end

% Compute relative odds ratios.
OR10=Wb10./Wb10(end,:);
OR5 = Wb5./Wb5(end,:);
OR1 = Wb1./Wb1(end,:);
OR0 = Wb0./Wb0(end,:);




% Plot.

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#852ED1','#AC70E0','#DB8624','#EFBD38','#2CD38F','#80E5B7','#eab3fa'};
GREY=[0.85,0.85,0.85];
Rs=getMscale(M); Rs=Rs*(100/max(Rs));

% Define expected Mlrg values.
n=1:Nc;
Ml=log10(1:Nc)/mean(b)+m1;
ql=0.10; qm=0.50; qh=0.90;
Ml1=Ml-log10(n.*(1-ql.^(1./n)))/mean(b);
Ml2=Ml-log10(n.*(1-qm.^(1./n)))/mean(b);
Ml3=Ml-log10(n.*(1-qh.^(1./n)))/mean(b);

% Plot EW-test results.
figure(7); clf;
% MvN plot.
ax1=subplot(311);
scatter([0,n]+1,[m1,M],[min(Rs),Rs],'r','filled','HandleVisibility','off'); hold on;
plot([0,n]+1,[m1,cummax(M)] ,'-r','HandleVisibility','off');
plot([0,n]+1,[m1,Ml2] ,'-.r','HandleVisibility','off');
plot([0,n]+1,[m1,Ml1] ,':r','HandleVisibility','off');
plot([0,n]+1,[m1,Ml3] ,':r','HandleVisibility','off');
for i=1:length(W10)
    plot((0:Nc)+1,[W10(i).Mmax(1),W10(i).Mmax],'-','Color',colours{i}, 'DisplayName',[name{i}]);
end
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 Nc+1]);
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvN plot.
ax2=subplot(312);
bp=bar([0,n]+1.5,Wb10,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(W10)
    plot([1 Nc+1],[i i]/length(W10),'-w');
    bp(i).FaceColor=colours{i};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 Nc+1]); ylim([0 1]);
set(gca, 'XScale', 'log');
% ORvN plot.
ax3=subplot(313);
for i=1:length(W10)-1
    loglog([0,n]+1,OR10(i,:),'-','Color',colours{i}, 'DisplayName',['M_{MAX} ',name{i}]); hold on;
end
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc+1]); ylim([1e-5 1e+5]);
linkaxes([ax1 ax2 ax3],'x');

% Plot order statistic depth comparison.
figure(59); clf;
% Zoomed out ORs.
subplot(211);
loglog([0,n]+1,max(OR10),'-b','LineWidth',1, 'DisplayName','M_{LRG(10)}'); hold on;
loglog([0,n]+1,max(OR5),'-.b','LineWidth',1, 'DisplayName','M_{LRG(5)}');
loglog([0,n]+1,max(OR1),'--b','LineWidth',1, 'DisplayName','M_{LRG(1)}');
loglog([0,n]+1,max(OR0), ':b','LineWidth',1, 'DisplayName','GR-MFD');
loglog(xlim(),1e0*[1 1],'-k', 'HandleVisibility','off');
loglog(xlim(),3e0*[1 1],':k', 'HandleVisibility','off');
loglog(xlim(),1e1*[1 1],':k', 'HandleVisibility','off');
loglog(xlim(),1e2*[1 1],':k', 'HandleVisibility','off');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc+1]); ylim([0.9 1.1*max(OR10(:))]);
legend('Location','northwest');
% Zoomed in ORs.
subplot(212);
loglog([0,n]+1,max(OR10),'-b','LineWidth',1, 'DisplayName','M_{LRG(10)}'); hold on;
loglog([0,n]+1,max(OR5),'-.b','LineWidth',1, 'DisplayName','M_{LRG(5)}');
loglog([0,n]+1,max(OR1),'--b','LineWidth',1, 'DisplayName','M_{LRG(1)}');
loglog([0,n]+1,max(OR0), ':b','LineWidth',1, 'DisplayName','GR-MFD');
loglog(xlim(),1e0*[1 1],'-k', 'HandleVisibility','off');
loglog(xlim(),3e0*[1 1],':k', 'HandleVisibility','off');
loglog(xlim(),1e1*[1 1],':k', 'HandleVisibility','off');
loglog(xlim(),1e2*[1 1],':k', 'HandleVisibility','off');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc+1]); ylim([1 1.1*1e3]);
legend('Location','northwest');








%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end