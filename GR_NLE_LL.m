function [LL,AICc,BIC]=GR_NLE_LL(Ma,Ms,m1,m2,a,b,i,T_flag,Wo)
  % Function that computes the log-likelihood (and similar scores) for the 
  % next i-th largest event, given a sample catalogue.
  % Code is not vectorized.  The numerical convolution integral broke this.
  % 
  % References:
  % Marzocchi & Sandri (2003). A review and new insights on the estimation of the b-valueand its uncertainty. Annals of geophysics.
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the upper magnitude bounds and length.
  n=length(Ma);
  K=3+1; % b,m1,m2.
  
  % Get the GR-MFD log-likelihood [Aki, 1965; Marzocchi & Sandri, 2003].
  [dm,pdf,~,~]=GR_NLE(Ma, m1,m2, a,b, i,T_flag,Wo);
  PDF=interp1(dm,pdf,Ms,'linear');
  LL=sum(log(PDF));
  
  % Compute the AIC & BIC statistics [Wagenmakers & Farrell, 2004].
  %AIC=2*K-2*LL;
  AICc=2*K+(2*K*(K+1)/(n-K-1))-2*LL;
  BIC=K*log(n)-2*LL;
  
end