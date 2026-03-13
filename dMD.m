function [dM,dMpdf,dMcdf,dMsvf]=dMD(pdf,d)
  % Function that computes the probability density function (PDF), 
  % cumulative distribution function (CDF), and survival function (SVF) of 
  % sequential variable differences dM - given the underlying PDF of M.
  % 
  % References:
  % Pishro-Nik (2014). Introduction to probability, statistics, and random processes. ISBN 0990637204.
  %
  % Written by Ryan Schultz.
  %
  
  % Estimate the distribution of dM, via a (numerical) convolution integral.
  [dMpdf,dM]=xcorr(pdf,'none'); % A cross-correlation is equivalent for dM=M2-M1.
  
  % Normalize and error handle the dM PDF.
  dMpdf=dMpdf*d*length(dMpdf)/length(dM(dM>=0));
  dMpdf=abs(dMpdf);
  
  % Clean up to the right units and domain.
  dM=dM*d;
  dMpdf(dM<0)=[]; dM(dM<0)=[];
  
  % Numerically estimate the CDF and SVF from the PDF.
  [~,dMcdf,dMsvf]=getDIST(dMpdf,d);
  
end




%%%% SUBROUNTINES.

% Given a PDF, find the CDF/SVF.
function [pdf,cdf,svf]=getDIST(pdf,dx)
    cdf=cumsum(pdf)*dx; %pdf=pdf/max(cdf);
    cdf=cdf/max(cdf);
    svf=1-cdf;
end