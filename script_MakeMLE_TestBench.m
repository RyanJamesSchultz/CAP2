% Test bench to explore the improved MLE-fitting.
% Used to make the data for Figure 6.
clear;

% Define some variables.
Ns=5e2;
Nc=[1e2 1e3];
dMexp=[-1.0 -0.8 -0.6 -0.4 -0.2 0.0 +1.0];
Mlrg=4.0;
b=1.0;
kd=10;
Nr=1e2;
type_flag='comp';

% Load in lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');

% Preallocate.
S=struct('Nc',0,'dMexp',0,'Mmax_bis',zeros([1 Ns]),'Mmax_std',zeros([1 Ns]));

% Loop over parameter space.
for j=1:length(Nc)
    for k=1:length(dMexp)
        
        % Save some important details.
        S(j,k).Nc=Nc(j);
        S(j,k).dMexp=dMexp(k);
        
        % Compute some dependent values.
        Mmax=Mlrg+dMexp(k);
        m1=Mlrg-log10(Nc(j))/b;
        m2=Mmax;
        a=log10(Nc(j));
        
        % Loop over all of the catalogue trials.
        for i=1:Ns
            [i,k,j]
            
            % Get the catalogue.
            M=GR_MFD_Rand(m1,m2, a,b, [1 Nc(j)]);
            
            % Fit m2.
            [m2_fit,nll]=M2fit(M,m1,b,kd,Nr,type_flag, Sa,ma,ba,m2a);
            m2_avg=mean(m2_fit);
            m2_std=std(m2_fit);
            
            % Stuff the results in the output data structure.
            S(j,k).Mmax_bis(i)=m2_avg-Mmax;
            S(j,k).Mmax_std(i)=m2_std;
        end
    end
end

% Save the data.
save('TestBenchMLE_temp.mat','S');


