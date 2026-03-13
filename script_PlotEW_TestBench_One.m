% Plots the test bench EW-test results, one example at a time.
% File upload size restriction prevented the upload of TestBenchEW_*.mat results, new files will need to be generated.
% See also script_MakeEW_TestBench.m.
clear;

% Load in lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');

% Load in the data from the test bench.
load('TestBenchEW_C10vC1.mat','S','name');
Ns=size(S,3);

% Pick which example(s) to examine further.
% i is the Nc index.
% j is the dMexp index.
% k is the catalogue realization index.
i=2; j=2; k=1:Ns;

% Get some relevant values.
Nc=S(i,j,1).Nc;
dMexp=S(i,j,1).dMexp;
type=S(i,j,1).type;
m1=S(i,j,1).m1;
b=S(i,j,1).b;

% Report.
type
dMexp
Nc

% Loop over all of the requested examples.
for ik=1:length(k)
    
    % Get the catalogue
    M=S(i,j,k(ik)).M;
    
    % Report KS-test details.
    k(ik)
    1-KS_dM_test(M,m1,'resample')
    %m1t=max([max(M)-2,m1]);
    %1-KS_dM_test(M,m1t,'resample')
    
    % Report MLE-test details.
    m2_fitA1=M2fit(M,m1,b,1,1e2,'approx');
    m2_fitB10=M2fit(M,m1,b,10,1e2,'comp',Sa,ma,ba,m2a);
    [mean(m2_fitA1)  std(m2_fitA1) ]
    [mean(m2_fitB10) std(m2_fitB10)]
    
    % Make weight matrices.
    W1=S(i,j,k(ik)).W1;
    W2=S(i,j,k(ik)).W2;
    Wb1=[]; Wb2=[];
    for l=1:length(W1)
        Wb1=[Wb1;W1(l).W];
        Wb2=[Wb2;W2(l).W];
    end
    OR1=Wb1./Wb1(end,:);
    OR2=Wb2./Wb2(end,:);
    
    % Define expected Mlrg values.
    n=1:Nc;
    Ml=log10(1:Nc)/mean(b)+m1;
    ql=0.10; qm=0.50; qh=0.90;
    Ml1=Ml-log10(n.*(1-ql.^(1./n)))/mean(b);
    Ml2=Ml-log10(n.*(1-qm.^(1./n)))/mean(b);
    Ml3=Ml-log10(n.*(1-qh.^(1./n)))/mean(b);
    
    
    
    % Plot.

    % Define some colours I'd like to use.
    colours={'#345da7','#587aff','#852ED1','#AC70E0','#DB8624','#EFBD38','#2CD38F','#80E5B7','#eab3fa'};
    GREY=[0.85,0.85,0.85];
    Rs=getMscale(M); Rs=Rs*(100/max(Rs));

    % Plot EW-test #1 results.
    figure(1); clf;
    % MvN plot.
    ax1=subplot(311);
    scatter((0:Nc)+1,[m1,M],[min(Rs),Rs],'r','filled','HandleVisibility','off'); hold on;
    plot((0:Nc)+1,[m1,cummax(M)] ,'-r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml2] ,'-.r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml1] ,':r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml3] ,':r','HandleVisibility','off');
    for l=1:length(W1)
        plot((0:Nc)+1,[W1(l).Mmax(1),W1(l).Mmax],'-','Color',colours{l}, 'DisplayName',[name{l}]);
    end
    xlabel('Chronological Event Number'); ylabel('Magnitude');
    xlim([1 Nc+1]); ylim([m1 max(M)+0.5])
    set(gca, 'XScale', 'log');
    legend('Location','northwest');
    title('EW-test #1');
    % WvN plot.
    ax2=subplot(312);
    bp=bar((0:Nc)+1.5,Wb1,'stacked','LineStyle','none','BarWidth',1); hold on;
    for l=1:length(W1)
        plot([1 Nc],[l l]/length(W1),'-w');
        bp(l).FaceColor=colours{l};
    end
    xlabel('Chronological Event Number'); ylabel('Model Weights');
    xlim([1 Nc]); ylim([0 1]);
    set(gca, 'XScale', 'log');
    % ORvN plot.
    ax3=subplot(313);
    for l=1:length(W1)-1
        loglog((0:Nc)+1,OR1(l,:),'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
    end
    loglog(xlim(),1*[1 1],'-k');
    loglog(xlim(),3*[1 1],':k');
    loglog(xlim(),10*[1 1],':k');
    loglog(xlim(),100*[1 1],':k');
    xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
    xlim([1 Nc]); ylim([1e-5 1e+10]);
    linkaxes([ax1 ax2 ax3],'x');
    
    % Plot EW-test #2 results.
    figure(2); clf;
    % MvN plot.
    ax1=subplot(311);
    scatter((0:Nc)+1,[m1,M],[min(Rs),Rs],'r','filled','HandleVisibility','off'); hold on;
    plot((0:Nc)+1,[m1,cummax(M)] ,'-r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml2] ,'-.r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml1] ,':r','HandleVisibility','off');
    plot([0,n]+1,[m1,Ml3] ,':r','HandleVisibility','off');
    for l=1:length(W2)
        plot((0:Nc)+1,[W2(l).Mmax(1),W2(l).Mmax],'-','Color',colours{l}, 'DisplayName',[name{l}]);
    end
    xlabel('Chronological Event Number'); ylabel('Magnitude');
    xlim([1 Nc+1]); ylim([m1 max(M)+0.5])
    set(gca, 'XScale', 'log');
    legend('Location','northwest');
    title('EW-test #2');
    % WvN plot.
    ax2=subplot(312);
    bp=bar((0:Nc)+1.5,Wb2,'stacked','LineStyle','none','BarWidth',1); hold on;
    for l=1:length(W2)
        plot([1 Nc],[l l]/length(W2),'-w');
        bp(l).FaceColor=colours{l};
    end
    xlabel('Chronological Event Number'); ylabel('Model Weights');
    xlim([1 Nc]); ylim([0 1]);
    set(gca, 'XScale', 'log');
    % ORvN plot.
    ax3=subplot(313);
    for l=1:length(W2)-1
        loglog((0:Nc)+1,OR2(l,:),'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
    end
    loglog(xlim(),1*[1 1],'-k');
    loglog(xlim(),3*[1 1],':k');
    loglog(xlim(),10*[1 1],':k');
    loglog(xlim(),100*[1 1],':k');
    xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
    xlim([1 Nc]); ylim([1e-5 1e+10]);
    linkaxes([ax1 ax2 ax3],'x');
    
    % Pause.
    pause;
end







%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end


