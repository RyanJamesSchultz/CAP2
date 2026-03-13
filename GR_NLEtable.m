function [PDFl]=GR_NLEtable(dM,m1,m2,a,b,i, S,dMa,ba,m2a)
  % Function that approximates the log-likelihood (using a lookup table) 
  % for the next i-th largest event, given a sample catalogue.
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Error handling for out of bounds Mmax values (like infinity).
  m2(m2>max(m2a))=max(m2a);
  
  % Get the NLE log-likelihood via a lookup table interpolation.
  PDFl=S(i).GI(dM,b,m2);
  
end