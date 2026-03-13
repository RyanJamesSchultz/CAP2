% Plots the test bench EW-test results averaged over all trials of one test.
% File upload size restriction prevented the upload of TestBenchEW_*.mat results, new files will need to be generated.
% See also script_MakeEW_TestBench.m.
clear;

% Load in the data from the test bench.
load('TestBenchEW_C10vC1.mat','S','name');

% Pick which example to examine further.
% i is the Nc index.
% j is the dMexp index.
i=2; j=2;

% Get some relevant values.
Ns=size(S,3);
Nm=length(S(i,j,1).W1);
Nc=S(i,j,1).Nc;
dMexp=S(i,j,1).dMexp;
type=S(i,j,1).type;

% Report values.
type
dMexp
Nc

% Get the averaged EW-test result.
Wb1=zeros([Nm,Nc+1]); Wb2=Wb1;
for k=1:Ns
    % Make weight matrices.
    W1=S(i,j,k).W1;
    W2=S(i,j,k).W2;
    Wb1i=[]; Wb2i=[];
    for l=1:Nm
        Wb1i=[Wb1i;W1(l).W];
        Wb2i=[Wb2i;W2(l).W];

        Wb1i(isnan(Wb1i))=0;
        Wb2i(isnan(Wb2i))=0;
    end
    Wb1=Wb1+(Wb1i);
    Wb2=Wb2+(Wb2i);
    
    % Make individual odds ratios.
    OR1i=Wb1./Wb1(end,:);
    OR2i=Wb2./Wb2(end,:);
    for l=1:Nm
        Sor1(l).or(k,:)=OR1i(l,:);
        Sor2(l).or(k,:)=OR2i(l,:);
    end
    
end
Wb1=(Wb1); Wb2=(Wb2);
Wb1=Wb1./sum(Wb1,1);
Wb2=Wb2./sum(Wb2,1);

% Get the index of the winning model.
OR1=Wb1./Wb1(end,:);
OR2=Wb2./Wb2(end,:);
[~,Im1]=max(OR1(:,end));
[~,Im2]=max(OR2(:,end));




% Plot.

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#852ED1','#AC70E0','#DB8624','#EFBD38','#2CD38F','#80E5B7','#eab3fa'};
GREY=[0.85,0.85,0.85];
%Rs=getMscale(M); Rs=Rs*(100/max(Rs));

% Plot average EW-test #1 results.
figure(1); clf;
% WvN plot.
ax1=subplot(311);
bp=bar((0:Nc)+1.5,Wb1,'stacked','LineStyle','none','BarWidth',1); hold on;
for l=1:length(W1)
    plot([1 Nc],[l l]/length(W1),'-w');
    bp(l).FaceColor=colours{l};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 Nc]); ylim([0 1]);
set(gca, 'XScale', 'log');
title('Average EW-test #1');
% ORvN plot.
ax2=subplot(312);
for l=1:length(W1)-1
    loglog((0:Nc)+1,OR1(l,:),'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
end
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc]); ylim([1e-5 1e+10]);
linkaxes([ax1 ax2],'x');
% ORvN plot 2.
ax2=subplot(313);
%for l=1:length(Sor1)
%    loglog((0:Nc)+1,Sor1(l).or,'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
%    loglog((0:Nc)+1,median(Sor1(l).or),'-','Color',colours{l},'LineWidth',3, 'HandleVisibility','off');
%end
loglog((0:Nc)+1,Sor1(Im1).or,'-','Color',colours{Im1}, 'DisplayName',['M_{MAX} ',name{Im1}]); hold on;
loglog((0:Nc)+1,median(Sor1(Im1).or),'-','Color','k','LineWidth',2, 'HandleVisibility','off');
loglog((0:Nc)+1,mean(Sor1(Im1).or),':','Color','k','LineWidth',2, 'HandleVisibility','off');
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc]); ylim([1e-5 1e+10]);
linkaxes([ax1 ax2],'x');

% Plot average EW-test #2 results.
figure(2); clf;
% WvN plot.
ax1=subplot(311);
bp=bar((0:Nc)+1.5,Wb2,'stacked','LineStyle','none','BarWidth',1); hold on;
for l=1:length(W2)
    plot([1 Nc],[l l]/length(W2),'-w');
    bp(l).FaceColor=colours{l};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 Nc]); ylim([0 1]);
set(gca, 'XScale', 'log');
title('Average EW-test #2');
% ORvN plot.
ax2=subplot(312);
for l=1:length(W2)-1
    loglog((0:Nc)+1,OR2(l,:),'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
end
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc]); ylim([1e-5 1e+10]);
linkaxes([ax1 ax2],'x');
% ORvN plot 2.
ax2=subplot(313);
%for l=1:length(Sor2)
%    loglog((0:Nc)+1,Sor2(l).or,'-','Color',colours{l}, 'DisplayName',['M_{MAX} ',name{l}]); hold on;
%    loglog((0:Nc)+1,median(Sor2(l).or),'-','Color',colours{l},'LineWidth',3, 'HandleVisibility','off');
%end
loglog((0:Nc)+1,Sor2(Im2).or,'-','Color',colours{Im2}, 'DisplayName',['M_{MAX} ',name{Im2}]); hold on;
loglog((0:Nc)+1,median(Sor2(Im2).or),'-','Color','k','LineWidth',2, 'HandleVisibility','off');
loglog((0:Nc)+1,mean(Sor2(Im2).or),':','Color','k','LineWidth',2, 'HandleVisibility','off');
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Nc]); ylim([1e-5 1e+10]);
linkaxes([ax1 ax2],'x');








%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end


