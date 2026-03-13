function [m2,nLL]=M2fit(M,m1,b,kd,Nr,type_flag, Sa,dMa,ba,m2a)
  % Function that MLE fits the best Mmax value, using magnitude differences 
  % in the sequence of the i-th largest events.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Set MLE fitting hyper parameters.
  options = optimoptions(@fmincon,'StepTolerance',1e-4,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-3,'Display','none');
  
  % Preallocate.
  m2=zeros([1 Nr+1]);
  nLL=m2;
  
  % Make M a row vector, then get catalogue size and the largest magnitude.
  M=M(:)';
  Nm=length(M);
  mlrg=max(M);

  % Format the b-values appropriately.
  if(length(b)==length(M))
      b=b(:)';
  else
      b=b*ones(size(M));
  end
  
  % Loop over all reshuffling trials.
  for l=1:Nr+1
      
      % Get a newly reshuffled catalogue.
      if(l==1)
          Ml=M;
      else
          Il=randperm(Nm);
          Ml=M(Il);
      end
      
      % Get the Mlrg & dMlrg sequences.
      D=OrderStatisticsI(Ml,kd,'unique'); 
      for i=1:kd
          D(i).I(isnan(D(i).Mlrg))=[];
          D(i).Mlrg(isnan(D(i).Mlrg))=[];
          D(i).Mlrg=[m1,D(i).Mlrg];
          D(i).dMlrg=diff(D(i).Mlrg);
          D(i).Mlrg(end)=[];

          % Error handling.
          if(isempty(D(i).dMlrg))
              disp('Warning: catalogue is too small for the user-input order statistic depth.');
              continue;
          end
      end
      
      % Find the optimal Mmax value (and corresponding log-likelihood).
      if(strcmpi(type_flag,'comp'))
          nll = @(m)  NllB( D,Ml,m1,m,b,kd, Sa,dMa,ba,m2a);
      elseif(strcmpi(type_flag,'approx'))
          nll = @(m)  NllA(D,Ml,m1,m,b,kd);
      elseif(strcmpi(type_flag,'both'))
          nll = @(m)  NllC(D,Ml,m1,m,b,kd, Sa,dMa,ba,m2a);
      end
      m2(l)=fmincon(nll,mlrg+0.5,[],[],[],[],mlrg,10,[],options);
      nLL(l)=nll(m2(l))-nll(Inf);
  end
  
end 




%%%% SUBROUNTINES.

% The negative log-likelihood for dMlrg(i).
function [nll]=NllB(D,M,m1,m2,b,kd, Sa,dMa,ba,m2a)
  
  % Start with the negative log-likelihood for the GR-MFD events.
  nll=-log(GR_MFD(M,m1,m2,0,b,'norm'));
  
  % Loop until the kd-th deepest order statistic.
  for i=1:kd
      
      % Add the negative log-likelihood for dMlrg(j).
      I=D(i).I;
      nll(I)=-GR_NLEtable(D(i).dMlrg,0,m2-D(i).Mlrg,0,b(I),i,Sa,dMa,ba,m2a);
  end
  nll=sum(nll);
  
end


% The approximate negative log-likelihood for dMlrg(i).
function [nll]=NllA(D,M,m1,m2,b,kd)
  
  % Start with the negative log-likelihood for the GR-MFD events.
  nll=-GR_MFD_LL(M,m1,m2,0,b);
  
  % Loop until the kd-th deepest order statistic.
  for i=1:kd
      
      % Add the (approximate) negative log-likelihood for dMlrg(j).
      I=D(i).I;
      nll=nll+GR_MFD_LL(D(i).dMlrg,0,m2-D(i).Mlrg,0,i*b(I));
  end
  
end


% The negative log-likelihood for dMlrg(i), combining both approaches.
function [nll]=NllC(D,M,m1,m2,b,kd, Sa,dMa,ba,m2a)
  
  % Start with the negative log-likelihood for the GR-MFD events.
  nll=-log(GR_MFD(M,m1,m2,0,b,'norm'));
  
  % Loop until the kd-th deepest order statistic.
  for i=1:kd
      
      % Add the negative log-likelihood for dMlrg(j).
      I=D(i).I;
      nll(I)=(nll(I)+log(GR_MFD( D(i).dMlrg,0,m2-D(i).Mlrg,0,i*b(I),'norm')))/2; % Approx.
      nll(I)=nll(I)-GR_NLEtable(D(i).dMlrg,0,m2-D(i).Mlrg,0,b(I),i,Sa,dMa,ba,m2a)/2; % Average with the composite form.
  end
  nll=sum(nll);
  
end



