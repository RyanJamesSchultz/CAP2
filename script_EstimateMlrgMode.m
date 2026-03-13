% Script to examine the general expression for the modal value of Mlrg(i).
% Used to show the analytical/numerical expressions of Mlrg(i) are the same.
% See also Supplementary Text S1.
clear;

% Define some constants.
b=1.0;
m1=0.0;
m2=Inf;
Nc=100;
d=1e-5;
k=[1:10,20:10:100];

% Derive some more variables.
m=m1:d:10;
a=log10(Nc);

% Expected distribution of M.
[Mpdf,Mcdf,Msvf]=GR_MFD(m, m1,m2, a,b, 'normalized'); Mpdf(m>m2)=0; Mcdf(m>m2)=1; Msvf(m>m2)=0;

% Loop over all of the choices for k.
Mlrg_num=zeros(size(k));
for i=1:length(k)
    
    % Expected distributions of Mlrg(i).
    [MLpdf,MLcdf,MLsvf]=OSD(Mpdf,Mcdf,Msvf,k(i),Nc);
    
    % Estimate the modal value of the Mlrg PDF.
    [~,j]=max(MLpdf);
    Mlrg_num(i)=m(j);
end

% Compute the analytical expression for Mlrg.
Mlrg_ana=log10(Nc./(Nc-k+1))/b+m1;



% Plot.

% Compare modal Mlrg values.
figure(1); clf;
subplot(211);
plot(k,(Mlrg_num-m1),'-o','DisplayName','Numerical'); hold on;
plot(k,(Mlrg_ana-m1),'--o','DisplayName','Analytical');
xlabel('Order Statistic, k'); ylabel('M_{LRG}');
legend('Location','northwest');
subplot(212);
semilogy(k,abs((Mlrg_num-Mlrg_ana)./Mlrg_ana),'-o'); hold on;
xlabel('Order Statistic, k'); ylabel('M_{LRG} residual');

% Plot the pdf and modal Mlrg value.
if(length(k)==1)
    YL=[1e-4 1e2];
    figure(2); clf;
    plot(m-m1,MLpdf,'-k'); hold on;
    plot(m(j)-m1,MLpdf(j),'or','MarkerFaceColor','r');
    plot(m-m1,100*MLsvf,'-k','LineWidth',2);
    xlabel('M-m_1'); ylabel('PDF/1-CDF');
    set(gca,'YScale','log');
    ylim(YL);
    %xlim([0 4]);
    title('M_{LRG}');
end


