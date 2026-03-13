function [LL,AICc,BIC]=GR_OSD_LL(M,m1,m2,a,b,k,N)
  % Function that computes the log-likelihood (and similar scores) for the 
  % k-th order statistic of the Gutenberg-Richter magnitude-frequency 
  % distributions (GR-MFD), given a sample catalogue.
  % Code is vectorized.
  % 
  % References:
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the upper magnitude bounds and length.
  n=length(M);
  K=3+1; % b,m1,m2.
  
  % Get the Mlrg(i,N) log-likelihood.
  [PDF,~,~]=GR_OSD(M, m1,m2, a,b, k,N);
  LL=sum(log(PDF));
  
  % Compute the AIC & BIC statistics [Wagenmakers & Farrell, 2004].
  %AIC=2*K-2*LL;
  AICc=2*K+(2*K*(K+1)/(n-K-1))-2*LL;
  BIC=K*log(n)-2*LL;
  
end