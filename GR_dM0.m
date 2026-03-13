function [PDF,CDF,SvF]=GR_dM0(dM,dm2,a,b,norm_flag)
  % Function that computes the distribution of differences in earthquake 
  % magnitudes from the Gutenberg-Richter magnitude-frequency distribution 
  % (GR-MFD).
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Define some useful values.
  cd=b.*log(10);
  xm=10.^(-b.*dM);
  cn=1-10.^(-b.*dm2);
  cm2=10.^(-2*b.*dm2);
  
  % The PDF & CDF for a doubly bounded GR-MFD.
  PDF=(cd.*xm.*(1-cm2./(xm.^2)))./(cn.^2);   PDF(dM>dm2)=NaN; PDF(dM<0)=NaN;
  CDF=(1-xm.*(1+cm2./(xm.^2))+cm2)./(cn.^2); CDF(dM>dm2)=NaN; CDF(dM<0)=NaN;
  SvF=(1-CDF);      SvF(dM>dm2)=NaN; SvF(dM<0)=NaN;
  
  % Change the normalization, if flagged to.
  if(strcmpi(norm_flag,'counts'))
      PDF=PDF.*10.^a;
      CDF=CDF.*10.^a;
      SvF=SvF.*10.^a;
  end
  
end