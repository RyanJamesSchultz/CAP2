function [S]=OrderStatisticsI(M,kd,truncation_flag)
  % Function that computes the first kd deepest order statistics of some 
  % sample.  This function assumes that the input sample is already 
  % chronologically sorted.
  % See also OrderStatistic.m.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Preallocate the output data structure.
  S=struct('Mlrg',[],'I',[]);
  
  % Loop over all requested order statistics.
  for i=1:kd
      
      % Get the largest values in sample.
      Mlrgi=cummax(M);
      
      % and then remove them from the sample.
      [Mlrgu,Ii]=unique(Mlrgi); S(kd).I=Ii;
      M(Ii)=NaN; % This will produce NaN values when there isn't an order statistic yet.
      
      % Save results into output structure (depending on input flag).
      if(strcmpi(truncation_flag,'unique'))
          S(i).Mlrg=Mlrgu;
      else
          S(i).Mlrg=Mlrgi;
      end
      S(i).dMlrg=diff(S(i).Mlrg);
      S(i).I=Ii;
  end
  
end