function [LL,AICc,BIC]=GR_NLE_LLtable(dM,m1,m2,a,b,i, Sa,dMa,ba,m2a)
  % Function that approximates the log-likelihood (using a lookup table) 
  % for the next i-th largest event, given a sample catalogue.
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 

  % Error handling for out of bounds Mmax values (like infinity).
  m2(m2>max(m2a))=max(m2a);
  
  % Get the upper magnitude bounds and length.
  n=length(dM);
  K=3+1; % b,m1,m2.
  
  % Get the NLE log-likelihood via a lookup table interpolation.
  %PDF=interp3(ba,dMa,m2a,(S(i).D), b,dM,m2, 'linear');  % Slow...
  PDFl=Sa(i).GI(dM,b,m2);
  LL=sum((PDFl));
  
  % Compute the AIC & BIC statistics [Wagenmakers & Farrell, 2004].
  %AIC=2*K-2*LL;
  AICc=2*K+(2*K*(K+1)/(n-K-1))-2*LL;
  BIC=K*log(n)-2*LL;
  
end