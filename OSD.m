function [OSpdf,OScdf,OSsvf]=OSD(pdf,cdf,svf,k,N)
  % Function that computes the probability density function (PDF), 
  % cumulative distribution function (CDF), and survival function (SVF) of 
  % a chosen order statistic (k,N) - given the underlying PDF, CDF, and SVF.
  % 
  % References:
  % David H.A., Nagaraja H.N. (2005) Order statistics. Wiley Series in Probability and Statistics, doi: 10.1002/0471722162.
  %
  % Written by Ryan Schultz.
  %
  
  % Slightly different function calls, depending on N (to better handle errors).
  if(N<=25)
      
      % Compute the PDF of the k-th order statistic [David & Nagaraja, 2005].
      OSpdf=gamma(N+1)/(gamma(k)*gamma(N-k+1))*pdf.*(cdf.^(k-1)).*(svf.^(N-k));
      
      % Compute the CDF (and SVF) of the k-th order statistic [David & Nagaraja, 2005].
      OScdf=zeros(size(OSpdf));
      for i=k:N
          OScdf=OScdf+gamma(N+1)/(gamma(i+1)*gamma(N-i+1))*(cdf.^(i)).*(svf.^(N-i));
      end
      OSsvf=1-OScdf;
      
  else
      
      % Compute the PDF of the k-th order statistic [David & Nagaraja, 2005].
      OSpdf=nchoosek(N,k)*k*pdf.*(cdf.^(k-1)).*(svf.^(N-k));
      
      % Compute the CDF (and SVF) of the k-th order statistic [David & Nagaraja, 2005].
      OScdf=zeros(size(OSpdf));
      for i=k:N
          OScdf=OScdf+nchoosek(N,k)*(cdf.^(i)).*(svf.^(N-i));
      end
      OSsvf=1-OScdf;
      
  end
  
end