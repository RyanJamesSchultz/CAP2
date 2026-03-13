function [PDF,CDF,SvF]=GR_OSD(M,m1,m2,a,b,k,N)
  % Function that computes the distribution of the k-th order statistic for 
  % the Gutenberg-Richter magnitude-frequency distribution (GR-MFD).
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the Mlrg(i,N) distribution.
  [pdf,cdf,svf]=GR_MFD(M, m1,m2, a,b, 'normalized');
  [PDF,CDF,SvF]=OSD(pdf,cdf,svf,k,N);
  
end