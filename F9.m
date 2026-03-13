% Plots the test bench EW-test results.
% File upload size restriction prevented the upload of TestBenchEW_*.mat results, new files will need to be generated.
% See also script_MakeEW_TestBench.m.
% Used to make Figures 9 & S10.
clear;

% Load in the data from the test bench.
load('TestBenchEW_B10vA1.mat','S','name'); S1=S; %F9.
%load('TestBenchEW_C10vA1.mat','S','name'); S1=S; % FS10.

% Pick which examples to examine further.
% i is the Nc index.
% j is the dMexp index.
i=2; j=1:7;
Im=5;

% Get some relevant values.
Ns=size(S,3);
Nm=length(S(i,j(1),1).W1);
Nc=S(i,j(1),1).Nc;
dMexp=[S(i,j,1).dMexp];
type={S(i,j,1).type};

% Report values.
type
dMexp
Nc

% Loop over all of the test bench dMexp cases.
for ij=1:length(j)

    % Loop over all of the test bench trials.
    Wb1=zeros([Nm,Nc+1]); Wb2=Wb1; Wb3=Wb1; Wb4=Wb1;
    for k=1:Ns
        
        % Make weight matrices.
        W1=S1(i,j(ij),k).W1;
        W2=S1(i,j(ij),k).W2;
        %W3=S2(i,j(ij),k).W1;
        %W4=S2(i,j(ij),k).W2;
        Wb1i=[]; Wb2i=[]; Wb3i=[]; Wb4i=[];
        for l=1:Nm
            Wb1i=[Wb1i;W1(l).W];
            Wb2i=[Wb2i;W2(l).W];
            %Wb3i=[Wb3i;W3(l).W];
            %Wb4i=[Wb4i;W4(l).W];
            
            Wb1i(isnan(Wb1i))=0;
            Wb2i(isnan(Wb2i))=0;
            %Wb3i(isnan(Wb3i))=0;
            %Wb4i(isnan(Wb4i))=0;
        end
        Wb1=Wb1+(Wb1i);
        Wb2=Wb2+(Wb2i);
        %Wb3=Wb3+(Wb3i);
        %Wb4=Wb4+(Wb4i);
        
        % Make individual odds ratios.
        OR1i=Wb1./Wb1(end,:);
        OR2i=Wb2./Wb2(end,:);
        %OR3i=Wb3./Wb3(end,:);
        %OR4i=Wb4./Wb4(end,:);
        for l=1:Nm
            Sor1(j(ij),l).or(k,:)=OR1i(l,:);
            Sor2(j(ij),l).or(k,:)=OR2i(l,:);
            %Sor3(j(ij),l).or(k,:)=OR3i(l,:);
            %Sor4(j(ij),l).or(k,:)=OR4i(l,:);
        end
    end
end




% Plot.

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#852ED1','#AC70E0','#DB8624','#EFBD38','#2CD38F','#80E5B7','#eab3fa'};

% Plot the odds ratio of the winning model, as a function of dMexp.
figure(9); clf;
semilogy(dMexp,ones(size(dMexp)),'-k'); hold on;
for l=Im
    [ORl,ORm,ORh]=getPrctiles(Sor1,j,l);
    semilogy(dMexp,ORm,'-o','Color',colours{l}); hold on;
    semilogy(dMexp,ORl,':','Color',colours{l});
    semilogy(dMexp,ORh,':','Color',colours{l});
    [ORl,ORm,ORh]=getPrctiles(Sor2,j,l);
    semilogy(dMexp,ORm,'-o','Color',colours{l}); hold on;
    semilogy(dMexp,ORl,':','Color',colours{l});
    semilogy(dMexp,ORh,':','Color',colours{l});
    %[ORl,ORm,ORh]=getPrctiles(Sor3,j,l);
    %semilogy(dMexp,ORm,'-o','Color',colours{l}); hold on;
    %semilogy(dMexp,ORl,':','Color',colours{l});
    %semilogy(dMexp,ORh,':','Color',colours{l});
    %[ORl,ORm,ORh]=getPrctiles(Sor4,j,l);
    %semilogy(dMexp,ORm,'-o','Color',colours{l}); hold on;
    %semilogy(dMexp,ORl,':','Color',colours{l});
    %semilogy(dMexp,ORh,':','Color',colours{l});
end
semilogy(xlim(),3*[1 1],':k');
semilogy(xlim(),10*[1 1],':k');
semilogy(xlim(),100*[1 1],':k');
xlabel('Degree-of-Truncation, \deltaM'); ylabel('Relative Odds Ratio');
ylim([1e-1 1e5]);






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Pl,Pm,Ph]=getPrctiles(S,j,Im)
  for ij=1:length(j)
      Pl(ij)=prctile(S(j(ij),Im).or(:,end),10);
      Pm(ij)=prctile(S(j(ij),Im).or(:,end),50);
      %Pm(ij)=mean(S(j(ij),Im).or(:,end));
      Ph(ij)=prctile(S(j(ij),Im).or(:,end),90);
  end
end
