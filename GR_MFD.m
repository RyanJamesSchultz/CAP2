function [PDF,CDF,SvF]=GR_MFD(M,m1,m2,a,b,norm_flag)
  % Function that computes the Gutenberg-Richter magnitude-frequency 
  % distributions (GR-MFD).
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Define some useful values.
  cd=b.*log(10);
  xm=10.^(-b.*(M-m1));
  cn=1-10.^(-b.*(m2-m1));
  
  % The PDF & CDF for a doubly bounded GR-MFD.
  PDF=(cd.*xm)./cn; PDF(M>m2)=NaN; PDF(M<m1)=NaN;
  CDF=(1  -xm)./cn; CDF(M>m2)=NaN; CDF(M<m1)=NaN;
  SvF=(1-CDF);      SvF(M>m2)=NaN; SvF(M<m1)=NaN;
  
  % Change the normalization, if flagged to.
  if(strcmpi(norm_flag,'counts'))
      PDF=PDF.*10.^a;
      CDF=CDF.*10.^a;
      SvF=SvF.*10.^a;
  end
  
end