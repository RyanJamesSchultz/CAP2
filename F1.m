% Script to show order statistics applied to a synthetic MvT sequence.
% Used to make Figure 1.
clear;

% Define some constants.
Nc=1e4;
Ns=1e2;
m1=0.0;
m2=3.5;
b=1.0;

% Derive some more variables.
T=1:Nc;
Mmax=linspace(0.8*(m2-m1)+m1,m2,Nc);
M=GR_MFD_Rand(m1,Mmax, log10(Nc),b, [1 Nc]);
Rs=getMscale(M)/3;

% Get the sequence of largest events (and deeper order statistics).
Mlrg1=OrderStatistic(M,length(M)-0,'none'); %Mlrg1=[m1,Mlrg1];
Mlrg2=OrderStatistic(M,length(M)-1,'none'); %Mlrg2=[m1,Mlrg2];
Mlrg3=OrderStatistic(M,length(M)-2,'none'); %Mlrg3=[m1,Mlrg3];
Mlrg4=OrderStatistic(M,length(M)-3,'none'); %Mlrg4=[m1,Mlrg4];

% Plot.
figure(1); clf;
scatter(T,M,Rs,'r','filled'); hold on;
plot(T,Mmax,'-b');
plot(T,Mlrg1 ,'-r');
plot(T,Mlrg2,'-m');
plot(T,Mlrg3,'-g');
plot(T,Mlrg4,'-c');
%set(gca, 'XScale', 'log');
xlabel('Chronological Event Index, i'); ylabel('Magnitude');
xlim([-max(T)/20 max(T)]);






%%%% SUBROUNTINE.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

