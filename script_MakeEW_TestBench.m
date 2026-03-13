% Test bench to explore the improved EW-test.
% Used to make the data for Figure 8.
% See also script_PlotEW_* for other ways to explore these datasets.
clear;

% Define some constants.
Ns=2e2;
Nc=[1e2 1e3];
dMexp=[   -1.0,   -0.8,   -0.6,   -0.4,   -0.2,    0.0,   +1.0,   +2.0,  -1.5,  -1.0,  -2/3,  -1/3];
Trunc={'const','const','const','const','const','const','const','const','logl','logl','logl','logl'};
Mlrg=4.0;
b=1.0;
dm=0.05;

% EW-test #1 parameters.
kd1=10;
Nr1=0;
reshuffle_flag1='none';
smooth_flag1='none';
GRtype_flag1='M';
NLEtype_flag1='both';

% EW-test #2 parameters.
kd2=1;
Nr2=0;
reshuffle_flag2='none';
smooth_flag2='none';
GRtype_flag2='M';
NLEtype_flag2='approx';

% Load in lookup table.
load('NLEtable_small.mat','Sa','ma','ba','m2a');

% Prep for Mmax guesses.
Kgr=2+1;

% Preallocate.
S=struct('Nc',0,'dMexp',0,'m1',0,'b',0,'type','','W1',[],'W2',[],'Sm',[],'M',[]);

% Loop over parameter space.
for j=1:length(Nc)
    for k=1:length(dMexp)
        
        % Save some important details.
        S(j,k,1).Nc=Nc(j);
        S(j,k,1).dMexp=dMexp(k);
        S(j,k,1).type=Trunc{k};
        
        % Compute some dependent values for the true Mmax.
        Mmax=Mlrg+dMexp(k);
        m1=Mlrg-log10(Nc(j))/b;
        a=log10(Nc(j));
        
        % Save more details.
        S(j,k,1).m1=m1;
        S(j,k,1).b=b;
        
        % Make the true Mmax function.
        if(strcmpi(Trunc{k},'const'))
            m2=Mmax;
        elseif(strcmpi(Trunc{k},'line'))
            m2=linspace(m1+0.1,Mmax,Nc(j));
        elseif(strcmpi(Trunc{k},'logl'))
            c0=m1+0.1;
            c1=(Mmax-c0)/log10(Nc(j));
            m2=c1*log10(1:Nc(j))+c0;
        end
        
        % Loop over all of the catalogue trials.
        for i=1:Ns
            [i,k,j]
            
            % Get the catalogue.
            M=GR_MFD_Rand(m1,m2, a,b, [1 Nc(j)]);
            cM=cummax(M);
            n=length(unique(cM));
            
            % Fit m2.
            %[m2_fit,nll]=M2fit(M,m1,b,kh,Nr,NLEtype_flag, Sa,ma,ba,m2a);
            %m2_avg=mean(m2_fit);
            %m2_std=std(m2_fit);
            
            % Make the Mmax guesses.
            name{1}='Linear';
            Ffxn = @(c1)  OF(cM(2:end),LinearM(2:Nc(j),M(1)+dm,c1),dm);
            c1=fmincon(Ffxn,20*max(M),[],[],[],[],0,Inf);
            Sm(1).Mmax=LinearM(1:Nc(j),M(1)+dm,c1);
            Sm(1).K=Kgr+2;
            name{2}='Linear2';
            Sm(2).Mmax=Sm(1).Mmax+0.1;
            Sm(2).K=Kgr+2;
            name{3}='Log-Linear';
            Ffxn = @(c1)  OF(cM(2:end),LoglinearM(2:Nc(j),M(1)+dm,c1),0);
            c1=fmincon(Ffxn,2*max(M),[],[],[],[],0,Inf);
            Sm(3).Mmax=LoglinearM(1:Nc(j),M(1)+dm,c1)+dm;
            Sm(3).K=Kgr+2;
            name{4}='Log-Linear2';
            Sm(4).Mmax=Sm(3).Mmax+0.1;
            Sm(4).K=Kgr+2;
            name{5}='Constant';
            Sm(5).Mmax=max(M)*ones(size(M))+0.1;
            Sm(5).K=Kgr+1;
            name{6}='Constant2';
            Sm(6).Mmax=max(M)+0.2;
            Sm(6).K=Kgr+1;
            name{7}='Envelope)';
            Sm(7).Mmax=cummax(M)+0.1;
            Sm(7).K=Kgr+2*n;
            name{8}='Envelope2';
            Sm(8).Mmax=cummax(M)+0.2;
            Sm(8).K=Kgr+2*n;
            name{9}='Unbound';
            Sm(9).Mmax=Inf;
            Sm(9).K=Kgr+0;
            
            % Do the EW-tests.
            W1=EnsembleW(M,m1,Sm,b,kd1,Nr1,reshuffle_flag1,smooth_flag1, GRtype_flag1,NLEtype_flag1,Sa,ma,ba,m2a);
            W2=EnsembleW(M,m1,Sm,b,kd2,Nr2,reshuffle_flag2,smooth_flag2, GRtype_flag2,NLEtype_flag2,Sa,ma,ba,m2a);
            
            % Stuff the results in the output data structure.
            S(j,k,i).W1=W1;
            S(j,k,i).W2=W2;
            S(j,k,i).Sm=Sm;
            S(j,k,i).M=M;
        end
    end
end

% Save the data.
save('TestBenchEW_temp.mat','S','name');













%%%% SUBROUNTINES.

% Linear spacing of M.
function [M]=LinearM(x,c0,c1)
  M=c1*x+(c0-c1);
end

% Log-Linear spacing of M.
function [M]=LoglinearM(x,c0,c1)
  M=c1*log10(x)+c0;
end

% Objective function for fitting Mlrg envelopes.
function [v]=OF(d,f,e)
  v=f-(d+e);
  if(any(v<0))
      df=f(2)-f(1);
      v=10*exp(-df);
  else
      v=min(v);
  end
end
