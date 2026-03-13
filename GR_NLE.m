function [dm,PDF,CDF,SvF]=GR_NLE(M,m1,m2,a,b,i,T_flag,Wo)
  % Function that computes the distribution of the next i-th largest event (NLE).
  % Code is not vectorized.  The numerical convolution integral breaks this.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Define the weight vector Wj.
  % See also script_DeepSumFactors.m.
  w=[];
  if(i==1)
      w=[1 0];
  elseif(i==2)
      w=[2.22 0.95 1];
  elseif(i==3)
      w=[1.75 1.51 0 1];
  elseif(i==4)
      w=[1.61 2.01 0 0 1];
  elseif(i==5)
      w=[4.94 0.00 5.66 0 0 1];
  elseif(i==6)
      w=[3.40 1.11 3.47 0 0 0 1];
  elseif(i==7)
      w=[5.50 0.00 7.83 0 0 0 0 1];
  elseif(i==8)
      w=[5.79 6.57 0.00 5.90 0 0 0 0 1];
  elseif(i==9)
      w=[10.93 8.34 0.0 13.08 0 0 0 0 0 1];
  elseif(i>=10)
      w=[10.14 0.00 10.38 5.79 0 0 0 0 0 0 1];
  end
  w=w/sum(w);
  
  % Define the fudgefactor.
  ff=[[0.00 0.42 0.65 0.81 0.92 1.01 1.08 1.14 1.20 1.24 1.28 1.33 1.34 1.38 1.40 1.42 1.46 1.48 1.50 1.50], 1.5*ones([1 20])];
  ff=ff(i);
  
  % Use input weights/fudgefactor, for fitting.
  if(Wo~=-1)
      if(strcmpi(T_flag,'fudge'))
          ff=Wo;
      else
          w=Wo;
          w=w/sum(w);
      end
  end
  
  % Use either the weights or fudgefactors.
  if(strcmpi(T_flag,'fudge'))
      
      % Get the 'fudged' distribution of the NLE
      [dm,PDF,CDF,SvF]=GR_dMD(M,m1,m2,a,b,1+ff,i+ff);
  else
      
      % Get the PDFs for the first possible order statistic.
      [dm,PDF,CDF,SvF]=GR_dMD(M,m1,m2,a,b,1,i);
      PDF=PDF*w(1); CDF=CDF*w(1); SvF=SvF*w(1);
      
      % Loop over deeper order statistics.
      for j=2:length(w)
          
          % Skip if this entry has zero weight.
          if(w(j)==0)
              continue;
          end
          
          % Get the PDFs for the j-th possible order statistic.
          [~,pdf,cdf,svf]=GR_dMD(M,m1,m2,a,b,j,i+j-1);
          PDF=PDF+pdf*w(j); CDF=CDF+cdf*w(j); SvF=SvF+svf*w(j);
      end
      
  end
  
end