% Script that determines when bounded/unbounded order statistics become indistinguishable.
% Used to inform Section 3.1 and make Figure 4.
clear;

% Define some constants.
Nc=1e4;
Ns=1e3;
m1=0.0;
b=1.0;
kd=10;
d=0.01;
type='norm';
O=-1;

% Derive some more variables.
m2a1=linspace(2.5,2.5,Nc); % Scenario A1.
m2a2=linspace(3.0,3.0,Nc); % Scenario A2.
m2a3=linspace(3.5,3.5,Nc); % Scenario A3.
m2b1=linspace(2.5,4.0,Nc); % Scenario B1.
m2b2=linspace(3.0,4.0,Nc); % Scenario B2.
m2b3=linspace(3.5,4.0,Nc); % Scenario B3.
a=log10(Nc);
mr=m1:d:10;

% Preallocate.
KSu=zeros([kd,Ns]);
KSbA1=KSu; MWxA1=KSu;
KSbA2=KSu; MWxA2=KSu;
KSbA3=KSu; MWxA3=KSu;
KSbB1=KSu; MWxB1=KSu;
KSbB2=KSu; MWxB2=KSu;
KSbB3=KSu; MWxB3=KSu;

% Loop over trials.
for i=1:Ns
    i
    
    % Draw a catalogue.
    Mu=GR_MFD_Rand(m1,Inf, a,b, [1 Nc]);
    MbA1=GR_MFD_Rand(m1,m2a1 , a,b, [1 Nc]);
    MbA2=GR_MFD_Rand(m1,m2a2 , a,b, [1 Nc]);
    MbA3=GR_MFD_Rand(m1,m2a3 , a,b, [1 Nc]);
    MbB1=GR_MFD_Rand(m1,m2b1 , a,b, [1 Nc]);
    MbB2=GR_MFD_Rand(m1,m2b2 , a,b, [1 Nc]);
    MbB3=GR_MFD_Rand(m1,m2b3 , a,b, [1 Nc]);
    
    % Loop until the deepest order statistic requested.
    for j=1:kd
        
        % Get the dMlrg(j) seqeunce (for the bounded catalogues).
        MlrgBa1=OrderStatistic(MbA1,length(MbA1)-(j-1), 'unique'); MlrgBa1=[m1,MlrgBa1]; MlrgBa1(isnan(MlrgBa1))=[];
        MlrgBa2=OrderStatistic(MbA2,length(MbA2)-(j-1), 'unique'); MlrgBa2=[m1,MlrgBa2]; MlrgBa2(isnan(MlrgBa2))=[];
        MlrgBa3=OrderStatistic(MbA3,length(MbA3)-(j-1), 'unique'); MlrgBa3=[m1,MlrgBa3]; MlrgBa3(isnan(MlrgBa3))=[];
        MlrgBb1=OrderStatistic(MbB1,length(MbB1)-(j-1), 'unique'); MlrgBb1=[m1,MlrgBb1]; MlrgBb1(isnan(MlrgBb1))=[];
        MlrgBb2=OrderStatistic(MbB2,length(MbB2)-(j-1), 'unique'); MlrgBb2=[m1,MlrgBb2]; MlrgBb2(isnan(MlrgBb2))=[];
        MlrgBb3=OrderStatistic(MbB3,length(MbB3)-(j-1), 'unique'); MlrgBb3=[m1,MlrgBb3]; MlrgBb3(isnan(MlrgBb3))=[];
        dM_lrgBa1=diff(MlrgBa1); dM_lrgBa2=diff(MlrgBa2); dM_lrgBa3=diff(MlrgBa3);
        dM_lrgBb1=diff(MlrgBb1); dM_lrgBb2=diff(MlrgBb2); dM_lrgBb3=diff(MlrgBb3);
        
        % Get the dMlrg(j) seqeunce (for the unbounded catalogue).
        MlrgU=OrderStatistic(Mu,length(Mu)-(j-1), 'unique'); MlrgU=[m1,MlrgU]; MlrgU(isnan(MlrgU))=[];
        dM_lrgU=diff(MlrgU);
        
        % Get the expected distributions of dMlrg(j).
        [dmU,dMLUpdf,dMLUcdf,dMLUsvf]=GR_NLE(mr, m1,    Inf, a,b, j,type,O);
        %[dmB,dMLBpdf,dMLBcdf,dMLBsvf]=GR_NLE(mr, m1,max(m2a1), a,b, j,type,O);
        
        % Draw a jump catalogue.
        [~,I]=unique(dMLUsvf);
        dMLu=interp1(dMLUsvf(I),dmU(I),rand([1 Nc]),'linear');
        
        % Test the sample distributions for similarity.
        [~,Pu ]=kstest2(dM_lrgU,dMLu); KSu(j,i)=log10(Pu);
        [~,PbA1 ]=kstest2(dM_lrgBa1,dMLu); KSbA1(j,i)=log10(PbA1);
        [~,PbA2 ]=kstest2(dM_lrgBa2,dMLu); KSbA2(j,i)=log10(PbA2);
        [~,PbA3 ]=kstest2(dM_lrgBa3,dMLu); KSbA3(j,i)=log10(PbA3);
        [~,PbB1 ]=kstest2(dM_lrgBb1,dMLu); KSbB1(j,i)=log10(PbB1);
        [~,PbB2 ]=kstest2(dM_lrgBb2,dMLu); KSbB2(j,i)=log10(PbB2);
        [~,PbB3 ]=kstest2(dM_lrgBb3,dMLu); KSbB3(j,i)=log10(PbB3);
        
    end
    
end

% Loops over the order statistics & trials.
for i=1:kd
    for j=1:Ns
        
        % Perform the Mann-Whitney U-test.
        I=randperm(Ns); I=I(1:round(0.9*length(I)));
        PxA1=ranksum(KSu(i,I),KSbA1(i,I),'tail','right'); MWxA1(i,j)=log10(PxA1);
        PxA2=ranksum(KSu(i,I),KSbA2(i,I),'tail','right'); MWxA2(i,j)=log10(PxA2);
        PxA3=ranksum(KSu(i,I),KSbA3(i,I),'tail','right'); MWxA3(i,j)=log10(PxA3);
        PxB1=ranksum(KSu(i,I),KSbB1(i,I),'tail','right'); MWxB1(i,j)=log10(PxB1);
        PxB2=ranksum(KSu(i,I),KSbB2(i,I),'tail','right'); MWxB2(i,j)=log10(PxB2);
        PxB3=ranksum(KSu(i,I),KSbB3(i,I),'tail','right'); MWxB3(i,j)=log10(PxB3);
    end
end




% Plot.

% The KS-test results, for the deepest order statistic.
figure(5); clf;
histogram(KSu(kd,:),round(sqrt(Ns)), 'FaceColor','k'); hold on;
histogram(KSbA1(kd,:),round(sqrt(Ns)), 'FaceColor','r');
plot(log10(0.05)*[1 1],ylim(), ':k');
xlabel('log_{10}(p-value)'); ylabel('Counts');

% Improvement over baseline.
figure(4); clf;
% KS-test ratios.
subplot(211);
m=mean(KSbA1-KSu,2); s=std(KSbA1-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#073a8c','MarkerFaceColor','#073a8c','DisplayName','Scenario A1'); hold on;
plot(m+s,1:kd,':','Color','#073a8c','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#073a8c','HandleVisibility','off');
m=mean(KSbA2-KSu,2); s=std(KSbA2-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#1c5fc9','MarkerFaceColor','#1c5fc9','DisplayName','Scenario A2');
plot(m+s,1:kd,':','Color','#1c5fc9','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#1c5fc9','HandleVisibility','off');
m=mean(KSbA3-KSu,2); s=std(KSbA3-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#3f82eb','MarkerFaceColor','#3f82eb','DisplayName','Scenario A3');
plot(m+s,1:kd,':','Color','#3f82eb','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#3f82eb','HandleVisibility','off');
m=mean(KSbB1-KSu,2); s=std(KSbB1-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#3b056e','MarkerFaceColor','#3b056e','DisplayName','Scenario B1');
plot(m+s,1:kd,':','Color','#3b056e','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#3b056e','HandleVisibility','off');
m=mean(KSbB2-KSu,2); s=std(KSbB2-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#7c09e8','MarkerFaceColor','#7c09e8','DisplayName','Scenario B2');
plot(m+s,1:kd,':','Color','#7c09e8','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#7c09e8','HandleVisibility','off');
m=mean(KSbB3-KSu,2); s=std(KSbB3-KSu,0,2);
plot(m+0,1:kd,'-o','Color','#9c45ed','MarkerFaceColor','#9c45ed','DisplayName','Scenario B3');
plot(m+s,1:kd,':','Color','#9c45ed','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#9c45ed','HandleVisibility','off');
plot([0 0],ylim(), ':k','HandleVisibility','off'); 
xlabel('log_{10}(p-value ratios)'); ylabel('Order Statistic Depth, i');
set(gca, 'YDir','reverse');
ylim([1 kd]);
legend('Location','southwest');
% MW histograms.
subplot(212);
m=mean(MWxA1,2); s=2*std(MWxA1,0,2);
plot(m+0,1:kd,'-o','Color','#073a8c','MarkerFaceColor','#073a8c','DisplayName','Scenario A1'); hold on;
plot(m+s,1:kd,':','Color','#073a8c','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#073a8c','HandleVisibility','off');
m=mean(MWxA2,2); s=std(MWxA2,0,2);
plot(m+0,1:kd,'-o','Color','#1c5fc9','MarkerFaceColor','#1c5fc9','DisplayName','Scenario A2');
plot(m+s,1:kd,':','Color','#1c5fc9','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#1c5fc9','HandleVisibility','off');
m=mean(MWxA3,2); s=std(MWxA3,0,2);
plot(m+0,1:kd,'-o','Color','#3f82eb','MarkerFaceColor','#3f82eb','DisplayName','Scenario A3');
plot(m+s,1:kd,':','Color','#3f82eb','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#3f82eb','HandleVisibility','off');
m=mean(MWxB1,2); s=2*std(MWxB1,0,2);
plot(m+0,1:kd,'-o','Color','#3b056e','MarkerFaceColor','#3b056e','DisplayName','Scenario B1');
plot(m+s,1:kd,':','Color','#3b056e','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#3b056e','HandleVisibility','off');
m=mean(MWxB2,2); s=std(MWxB2,0,2);
plot(m+0,1:kd,'-o','Color','#7c09e8','MarkerFaceColor','#7c09e8','DisplayName','Scenario B2');
plot(m+s,1:kd,':','Color','#7c09e8','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#7c09e8','HandleVisibility','off');
m=mean(MWxB3,2); s=std(MWxB3,0,2);
plot(m+0,1:kd,'-o','Color','#9c45ed','MarkerFaceColor','#9c45ed','DisplayName','Scenario B3');
plot(m+s,1:kd,':','Color','#9c45ed','HandleVisibility','off');
plot(m-s,1:kd,':','Color','#9c45ed','HandleVisibility','off');
plot(log10(0.05)*[1 1],ylim(), '--k','HandleVisibility','off');
xlabel('log_{10}(p-value)'); ylabel('Order Statistic Depth, i');
set(gca, 'YDir','reverse');
ylim([1 kd]);
legend('Location','southwest');
