function [dm,PDF,CDF,SvF]=GR_dMD(M,m1,m2,a,b,k,N)
  % Function that numerically computes the distributions of k-th order 
  % statistic for magnitude differences in the Gutenberg-Richter 
  % magnitude-frequency distribution (GR-MFD).
  % Code is not vectorized.  The numerical convolution integral broke this.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the sampling difference.
  d=M(2)-M(1);

  % Get the dMlrg(i,N) distribution.
  [Mpdf,Mcdf,Msvf]=GR_MFD(M, m1,m2, a,b, 'normalized'); Mpdf(M>m2)=0; Mcdf(M>m2)=1; Msvf(M>m2)=0;
  [Opdf,Ocdf,Osvf]=OSD(Mpdf,Mcdf,Msvf,k,N);
  [dm,PDF,CDF,SvF]=dMD(Opdf,d);
  
end