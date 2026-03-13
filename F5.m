% Script that visualizes the negative log-likelihood for Mmax fitting.
% Used to make Figures 5 & S4-S7.
clear;
rng(2); % Regular case (must change m2 to 2.3; Figures 5, S4, S6, & S7).
%rng(8); % Unbound case (must change m2 to Inf; Figure S5).

% Define some constants.
Nc=1e3;
Nr=1e2;
m1=0.0;
b=1.0;
kd=10;
d=0.01;

% Derive some more variables.
m2=linspace(2.3,2.3,Nc);
a=log10(Nc);
Mmax=max(m2);

% Get a catalogue.
M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);
M2=max(M):d:max(M)+2;

% Load in lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');

% Preallocate.
D1=struct('Mlrg',[],'dMlrg',[],'nLL',[],'nLLA',[],'nLLB',[]);

% Loop over all catalogue reshuffle trials.
for l=1:Nr
    l
    
    % Get a newly reshuffled catalogue.
    if(l==1)
        In=1:Nc;
    else
        In=randperm(Nc);
    end
    Mi=M(In);
    
    % Collect the observable order statistics.
    D2=OrderStatisticsI(Mi,kd,'unique');
    for i=1:kd
        Mlrgi=D2(i).Mlrg;
        Mlrgi=[m1,Mlrgi]; Mlrgi(isnan(Mlrgi))=[]; dMlrgi=diff(Mlrgi); Mlrgi(end)=[];
        D1(i).Mlrg=Mlrgi;
        D1(i).dMlrg=dMlrgi;
        D1(i).I=D2(i).I;
        D1(i).I(isnan(D2(i).Mlrg))=[];
    end
    %m2_fitA(l)=M2fit(Mi,m1,b,kh,0,'approx');
    %m2_fitB(l)=M2fit(Mi,m1,b,kh,0,'comp', S,ma,ba,m2a);
    
    % Compute the differences in magnitudes.
    %dMi=abs(diff([m1,M]));

    % Loop over all of the possible Mmax values.
    for j=1:length(M2)
        
        % Compute the GR-MFD negative log-likelihood.
        D1(kd+1).nLL(j,l)=-GR_MFD_LL(Mi,m1,M2(j),a,b)+GR_MFD_LL(Mi,m1,Inf,a,b);
        %nll=-log(GR_MFD(Mi,m1,M2(j),a,b,'norm'))+log(GR_MFD(Mi,m1,Inf,a,b,'norm'));
        %log(GR_MFD(M,m1,m2,0,b,'norm'));

        % Compute the dM negative log-likelihood.
        %D1(kh+2).nLL(j,l)=-GR_dM0_LL(dMi,M2(j),a,b)+GR_dM0_LL(dMi,Inf,a,b);
        
        % Loop over all of the order statistics.
        for i=1:kd
            Mlrgi=D1(i).Mlrg; dMlrgi=D1(i).dMlrg;
            I=D1(i).I;
            
            % Compute the dMlrg(i) negative log-likelihood.
            D1(i).nLLA(j,l)=+GR_MFD_LL(dMlrgi,0,M2(j)-Mlrgi,a,i*b)-GR_MFD_LL(dMlrgi,0,Inf,a,i*b); % Approximation.
            %D1(i).nLLB(j,l)=-GR_NLE_LLtable(dMlrgi, 0,M2(j)-Mlrgi, a,b*ones(size(Mlrgi)), i, S,ma,ba,m2a) + GR_NLE_LLtable(dMlrgi, 0,Inf*ones(size(Mlrgi)), a,b*ones(size(Mlrgi)), i, S,ma,ba,m2a);
            D1(i).nLLB(j,l)=-GR_NLE_LLtable(dMlrgi, 0,M2(j)-Mlrgi, a,b*ones(size(Mlrgi)), i, Sa,ma,ba,m2a) + GR_NLE_LLtable(dMlrgi, 0,Inf*ones(size(Mlrgi)), a,b*ones(size(Mlrgi)), i, Sa,ma,ba,m2a)+GR_MFD_LL(Mi(I),m1,M2(j),a,b)-GR_MFD_LL(Mi(I),m1,Inf,a,b);

            %nll(I)=-GR_NLEtable(dMlrgi,0,M2(j)-Mlrgi,0,b*ones(size(Mlrgi)),i,S,ma,ba,m2a)+GR_NLEtable(dMlrgi,0,Inf*ones(size(Mlrgi)),0,b*ones(size(Mlrgi)),i,S,ma,ba,m2a);
            %D1(i).nLLB(j,l)=sum(nll);
        end
        %D1(i).nLLB(j,l)=sum(nll);
        %sum(nll)
    end
end
%%
% Preallocate.
N1=struct('nLL',[],'Mfit',[],'Mfit_avg',[],'nLL_min',[]);

% Loop over all the order statistics.
nLLgr=D1(kd+1).nLL;
nLLA=zeros(size(D1(kd+1).nLL));
nLLB=nLLA;
for i=1:kd
    
    % Get the cumulative negative log-likelihood function.
    nLLA=nLLA+D1(i).nLLA;
    nLLB=nLLB+D1(i).nLLB;
    
    % Decide which of the two nLL options will be used for plotting.
    %nLL=(2/2)*nLLA + (0/2)*nLLB + (2/2)*nLLgr; % Approximate.
    nLL=(0/2)*nLLA + (2/2)*nLLB + (2/2)*nLLgr; % Composite form.
    %nLL=(1/2)*nLLA + (1/2)*nLLB + (2/2)*nLLgr; % Average of the two.
    N1(i).nLL=nLL;
    
    % Get the MLE estimates of Mmax.
    [nLL_min,Im]=min(nLL,[],1);
    N1(i).Mfit=M2(Im);
    N1(i).Mfit_avg=mean(N1(i).Mfit);
    N1(i).nLL_min=mean(nLL_min);
end

% Get the minimum nLL values.
for i=1:length(Im)
    LLfit(i)=nLL(Im(i),i);
end

% Report values.
max(M)
Mmax
[mean(N1(kd).Mfit) median(N1(kd).Mfit)]
std(N1(kd).Mfit)




% Plot.
colours={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#A2142F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#A2142F'};

% nLL pieces for each order statistic.
figure(5); clf;
% Cumulative negative log-likelihoods.
subplot(211);
plot(M2,mean(D1(kd+1).nLL,2),'DisplayName','GR-MFD'); hold on;
for i=1:kd
    y=mean(N1(i).nLL,2);
    plot(M2,y,'Color',colours{i+1},'DisplayName',['+\DeltaM_{LRG(',num2str(i),')}']);
    plot(N1(i).Mfit_avg,interp1(M2,y,N1(i).Mfit_avg,'linear'),'o','MarkerFaceColor',colours{i+1},'MarkerEdgeColor','k','HandleVisibility','off');
end
plot(Mmax*[1 1],ylim(),'-k','HandleVisibility','off');
plot(xlim(),[0 0],':k','HandleVisibility','off');
xlabel('Possible M_{MAX} Estimates'); ylabel('Cumulative Negative Log-Likelihood');
xlim([min(M2) max(M2)]); 
Yl=ylim(); ylim([Yl(1) -Yl(1)]);
%ylim([-20 +20]);
legend('location','northeast');
% Plot the histogram of MLE fitted Mmax values.
subplot(212);
histogram(max(M)*ones([1 Nr]),M2); hold on;
for i=[1 kd]
    histogram(N1(i).Mfit,M2,'FaceColor',colours{i+1});
    plot(N1(i).Mfit_avg*[1 1],ylim(),'--','LineWidth',2,'Color',colours{i+1});
end
plot(Mmax*[1 1],ylim(),'-k');
xlabel('Possible M_{MAX} Estimates'); ylabel('Counts');
xlim([min(M2) max(M2)]); ylim([0 Nr/2]);

% Best fit to the cumulative nLL.
figure(54); clf;
% Plot the negative log-likelihood for each reshuffle.
subplot(211);
plot(M2,nLL,'-r'); hold on;
plot(M2,nLL(:,1),'-b');
plot(M2,mean(nLL,2),'-k','LineWidth',2);
plot(N1(kd).Mfit_avg*[1 1],ylim(),'--r','LineWidth',2);
plot(Mmax*[1 1],ylim(),'-k');
plot(xlim(),[0 0],':k');
xlabel('Possible M_{MAX} Estimates'); ylabel('Negative Log-Likelihood');
xlim([min(M2) max(M2)]); Yl=ylim(); ylim([Yl(1) -Yl(1)]);
% Plot the histogram of MLE fitted Mmax values.
subplot(212);
histogram(N1(kd).Mfit,M2,'FaceColor','r'); hold on;
plot(N1(kd).Mfit_avg*[1 1],ylim(),'--r','LineWidth',2);
plot(Mmax*[1 1],ylim(),'-k');
xlabel('Possible M_{MAX} Estimates'); ylabel('Counts');
xlim([min(M2) max(M2)]);

% Extra plot for building intuition.
figure(1); clf;
plot(M2,(mean(D1(kd+1).nLL,2)),'DisplayName','GR-MFD'); hold on;
for i=1:kd
    y=mean(D1(i).nLLB,2);
    plot(M2,(y),'--','Color',colours{i+1},'DisplayName',['Full \DeltaM_{LRG(',num2str(i),')}']);
    y=mean(D1(i).nLLA,2);
    plot(M2,(y),':','Color',colours{i+1},'DisplayName',['Approx \DeltaM_{LRG(',num2str(i),')}']);
end
xlabel('Possible M_{MAX} Estimates'); ylabel('Cumulative Negative Log-Likelihood');
xlim([min(M2) max(M2)]);
plot(xlim(), [0 0],'-k','HandleVisibility','off');
legend('location','northeast');

% See the GR-MFD.
figure(10); clf; 
semilogy(sort(M,'descend'),1:length(M),'-b'); hold on; 
x=min(M):0.01:max(M);
[~,~,svf1]=GR_MFD(x,m1,max(m2),a,b,'count');
[~,~,svf2]=GR_MFD(x,m1,Inf,a,b,'count');
semilogy(x,svf1,'-k');
semilogy(x,svf2,':k');



