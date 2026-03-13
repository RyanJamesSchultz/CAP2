% Script to apply EW-test to the Groningen data.
% Used to make Figure 10.
clear;

% Define some values.
Ns=5e1;
m1=1.5;
dM=0.1;
kd=1;
Rflag1='resample';
Rflag2='regular';
Nr1=1e2;
Nr2=1e2;
SMOOTHflag='mean';
GRtype_flag='M';
NLEtype_flag='both';
Kgr=2+1;
q=-1;

% Set the RNG seed in stone, so the plot is reproducible.
rng(1);

% Load in the NLE lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');

% Get the Groningen catalogue.
load('data/Groningen.mat','Catalog'); 
ID=[Catalog(1).val]'; Lat=[Catalog(3).val]'; Lon=[Catalog(4).val]'; Dep=[Catalog(5).val]';
T=[Catalog(2).val]'; T=datetime(T,'ConvertFrom','datenum');
M=[Catalog(6).val]';

% Fitler the catalogue spatiotemporally.
[ID,Lon,Lat,Dep,T,M,~]=filtCat(ID,Lon,Lat,Dep,T,M);

% Loop over the bootstrap trials.
for i=1:Ns
    i
    
    % Apply a perturbation.
    Mi=M+(dM*rand(size(M))-dM/2);
    m1i=m1+((dM/2)*rand(size(m1))-dM/4);
    
    % Truncate on the magnitude of completeness.
    Im=(Mi>=m1i);
    %ID=ID(Im); Lat=Lat(Im); Lon=Lon(Im); Dep=Dep(Im); T=T(Im); M=M(Im);
    
    % Derive some variables.
    Nc=length(Mi(Im));
    
    % Fit the GR-MFD b-value.
    [bi,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(Mi,m1i,dM);
    b(i)=bi;
    
    % Get the expected Mlrg value.
    Mlrg_est(i)=m1i+log10(Nc)/bi;
    
    % Do the MLE-test.
    [m2_fitA,nllA]=M2fit(Mi(Im),m1i,bi,kd, Nr1,'approx');
    [m2_fitB,nllB]=M2fit(Mi(Im),m1i,bi,kd,Nr1,'comp', Sa,ma,ba,m2a);
    [m2_fitC,nllC]=M2fit(Mi(Im),m1i,bi,kd,Nr1,'both', Sa,ma,ba,m2a);
    m2A_avg(i)=mean(m2_fitA);
    m2B_avg(i)=mean(m2_fitB);
    m2C_avg(i)=mean(m2_fitC);
    m2A_std(i)=std(m2_fitA);
    m2B_std(i)=std(m2_fitB);
    m2C_std(i)=std(m2_fitC);
    
    % Do the KS-test.
    Msamp=GR_MFD_Rand(m1i,m2A_avg(i), log10(Nc),bi, [1 Nc]);
    [~,p_comp]=kstest2(Mi(Im),Msamp);
    KSp_fl(i)=p_comp;
    KSp_dm(i)=KS_dM_test(Mi(Im),m1i,Rflag1);
    
end

% Truncate on the magnitude of completeness.
Im=(M>=m1);
Nc=length(M(Im));

% Get the sequence of largest events (and indicies to when they occur).
Mlrg=OrderStatistic(M(Im),Nc-0,'none');
[~,Id]=unique(Mlrg);
Tlrg=T(Im);

% Propose possible three possible models of Mmax.
Sm(1).Mmax=(mean([m2A_avg,m2B_avg,m2C_avg]))*ones(size(Mlrg));
Sm(1).K=Kgr+1;
Sm(2).Mmax=(mean([m2A_avg,m2B_avg,m2C_avg])+mean(m2A_std))*ones(size(Mlrg));
Sm(2).K=Kgr+1;
Sm(3).Mmax=(Inf)*ones(size(Mlrg));
Sm(3).K=Kgr+0;

% Do the EW-test.
W=EnsembleW(M(Im),m1,Sm,mean(b),kd,Nr2,Rflag2,SMOOTHflag, GRtype_flag,NLEtype_flag,Sa,ma,ba,m2a);

% Get each Mmax model's estimate of the NLE's magnitude.
MnleE=zeros(size(Mlrg)); Wb=[];
for j=1:length(W)
    W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,mean(bi),q);
    MnleE=MnleE+(W(j).Mnle.*W(j).W(2:end));
    Wb=[Wb;W(j).W];
end

% Get the odds ratios, relative to the unbound model.
OR=Wb./Wb(end,:);

% Get the GR-MFD plotting values.
po=[-mean(b),a];
Mgr_fit=[m1, max(M)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Report the simple-test results.
disp('Simple-tests');
[max(M) mean(Mlrg_est) max(M)-mean(Mlrg_est)]
sum(M>m1)

% Report the KS-test p-value.
disp('KS(ΔM)-tests');
geomean(KSp_dm)
geomean(KSp_fl)

% Report the MLE-test Mmax values.
disp('MLE(ΔM)-tests');
[mean(m2A_avg) mean(m2A_std)]
[mean(m2B_avg) mean(m2B_std)]

% Report the EW-test model weights and the relative odds. 
disp('EW(ΔM)-tests');
Wb(:,end)
[Wb(:,end)/Wb(end,end)]




% Plot.

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#eab3fa'};
styles={'-','--','-'};
names={'Mean','Mean + error','Unbound'};
GREY=[0.85,0.85,0.85];
Rs=getMscale(M)*200; Rs=Rs*(100/max(Rs));

% Define expected Mlrg values.
n=1:Nc;
Ml=log10(1:Nc)/mean(b)+m1;
ql=0.10; qm=0.50; qh=0.90;
Ml1=Ml-log10(n.*(1-ql.^(1./n)))/mean(b);
Ml2=Ml-log10(n.*(1-qm.^(1./n)))/mean(b);
Ml3=Ml-log10(n.*(1-qh.^(1./n)))/mean(b);

% Plot time series data.
figure(1); clf;
% MvT plot.
scatter(T,M,Rs,'r','filled','HandleVisibility','off'); hold on;
plot(T,cummax(M),'-r');
plot(xlim(),m1*[1 1],'--r');
xlabel('Time'); ylabel('Magnitude');
ylim([m1-0.2 max(M)]);

% Plot EW-test results.
figure(10); clf;
% MvT plot.
ax1=subplot(311);
scatter([0,n]+1,[m1,M(Im)],[min(Rs(Im)),Rs(Im)],'r','filled','HandleVisibility','off'); hold on;
plot([0,n]+1,[m1,Ml2] ,'-.r','HandleVisibility','off');
plot([0,n]+1,[m1,Ml1] ,':r','HandleVisibility','off');
plot([0,n]+1,[m1,Ml3] ,':r','HandleVisibility','off');
plot([0,n]+1,[m1,Mlrg],'-r','HandleVisibility','off');
for i=1:length(W)
    plot([0,n]+1,[W(i).Mmax(1),W(i).Mmax],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{MAX} ',names{i}]);
    %plot([0,n]+1,[W(i).Mnle(1),W(i).Mnle],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{NLE} ',names{i}]);
end
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 Nc+1]); ylim([m1, max(ylim())]);
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvT plot.
ax2=subplot(312);
bp=bar([0,n]+1.5,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(W)
    plot([1 Nc+1],[i i]/length(W),'-w');
    bp(i).FaceColor=colours{i};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1,Nc+1]); ylim([0 1]);
set(gca, 'XScale', 'log');
% ORvN plot.
ax3=subplot(313);
for i=1:length(W)-1
    loglog([0,n]+1,OR(i,:),'-','Color',colours{i}, 'DisplayName',['M_{MAX} ',names{i}]); hold on;
end
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc+1]); ylim([1e-1 1e+3]);
linkaxes([ax1 ax2 ax3],'x');

% Plot catalogue info.
figure(3); clf;
% Map.
ax1=subplot(3,3,[1 2 4 5]);
scatter(Lon(~Im),Lat(~Im),Rs(~Im),'m','filled'); hold on;
scatter(Lon(Im), Lat(Im), Rs(Im), 'r','filled');
xlabel('Longitude'); ylabel('Latitude');
% Depth.
ax2=subplot(3,3,[3 6]);
scatter(Dep(~Im),Lat(~Im),Rs(~Im),'m','filled'); hold on;
scatter(Dep(Im), Lat(Im), Rs(Im), 'r','filled');
xlabel('Depth (km)'); ylabel('Latitude');
% GR-FMD.
ax3=subplot(3,3,[7 8 9]);
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);
plot(m1*[1 1],ylim,'--k');
xlabel('Magnitude'); ylabel('Count');






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end


% Filter the catalogue spatiotemporally.
function [id,x,y,z,t,m,I,m1]=filtCat(id,x,y,z,t,m)
  %I=0;
  
  % Define filter details.
  z_L=[0 10];
  x_L=[ 6.5  7.0  7.0  6.5];
  y_L=[53.5 53.5 53.1 53.1];
  T_L=[datetime(1980,01,01) datetime(2023,01,01)];

  % Filter laterally, temporally and in depth.
  Is=inpolygon(x,y,x_L,y_L);
  It=(t>=min(T_L))&(t<=max(T_L));
  Id=(z>=min(z_L))&(z<=max(z_L));
  
  % Apply the filter.
  I=Is&Id&It;
  id=id(I);
  x=x(I);
  y=y(I);
  z=z(I);
  t=t(I);
  m=m(I);
end

