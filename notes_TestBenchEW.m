% Notes reagrding the differences between EW-test types, inferred from the
% test bench results.
% 
% Terminology:
% A - Approximate form of nLL.
% B - Composite form of nLL.
% C - Best/combined form of nLL.
% 
% General notes.
% 10-100s of events seems to be enough to distinguish bound from unbound.
% 100-1000 of events seems to be enough to discern the functional form of Mmax.
% 
% TestBenchEW_B10vA10.mat
% B10 & A10 parameters.
% kh1=10; Nr1=0; reshuffle_flag1='none'; smooth_flag1='mean'; GRtype_flag1='M'; NLEtype_flag1='comp';
% kh2=10; Nr2=0; reshuffle_flag2='none'; smooth_flag2='mean'; GRtype_flag2='M'; NLEtype_flag2='approx';
% B10 is better at discerning bound cases.
% B10 can discern constant up to dMexp=0.  Wow!
% A10 is better at discerning unbound cases.
% B10 is more prone to sometimes getting (individual) unbound cases wrong.
% B10 is better at discerning the correct functional form of Mmax.
% B10 can sometimes initially pick the wrong functional form of Mmax.
% 
% TestBenchEW_B10vC10.mat
% B10 & C10 parameters.
% kh1=10; Nr1=0; reshuffle_flag1='none'; smooth_flag1='mean'; GRtype_flag1='M'; NLEtype_flag1='comp';
% kh2=10; Nr2=0; reshuffle_flag2='none'; smooth_flag2='mean'; GRtype_flag2='M'; NLEtype_flag2='both';
% C10 generally takes more samples to make an initial assessment.
% C10 generally becomes similarly confident of final assessment, on average. 
% C10 is better at identifying the true functional form of Mmax.
% C10 is better at discerning unbound cases.
% C10 seems to be a decent compromise, likely a good starting point for any given catalogue.
% 
% TestBenchEW_C10vC10rf.mat
% C10 & C10rf parameters.
% kh1=10; Nr1=0;  reshuffle_flag1='none';  smooth_flag1='mean'; GRtype_flag1='M'; NLEtype_flag1='both';
% kh2=10; Nr2=10; reshuffle_flag2='fixed'; smooth_flag2='mean'; GRtype_flag2='M'; NLEtype_flag2='both';
% Fixed reshuffling has limited impact on asessing unbound cases.
% Fixed reshuffling generally improves ability to discern bound cases.
% Fixed reshuffling generally imporves ability to discern the correct functional form of Mmax.
% 
% TestBenchEW_C10vC10rb.mat
% C10 & C10rb parameters.
% kh1=10; Nr1=0;  reshuffle_flag1='none';    smooth_flag1='mean'; GRtype_flag1='M'; NLEtype_flag1='both';
% kh2=10; Nr2=10; reshuffle_flag2='blocked'; smooth_flag2='mean'; GRtype_flag2='M'; NLEtype_flag2='both';
% Significantly improves confidence for bound cases.
% Exacerbates B10's problems in incorrectly discerning unbound cases.
%
% TestBenchEW_C10vC10rr.mat
% C10 & C10rr parameters.
% kh1=10; Nr1=0;  reshuffle_flag1='none';    smooth_flag1='mean'; GRtype_flag1='M'; NLEtype_flag1='both';
% kh2=10; Nr2=10; reshuffle_flag2='regular'; smooth_flag2='mean'; GRtype_flag2='M'; NLEtype_flag2='both';
% Doesn't generally work for Mmax(t) approaches.  Makes impossible scenarios.
% Drastically improves ability to discern (constant) bound cases.
% 
% TestBenchEW_B10vB10dm.mat
% B10 & B10dm parameters.
% kh1=10; Nr1=0; reshuffle_flag1='none'; smooth_flag1='mean'; GRtype_flag1='M';  NLEtype_flag1='comp';
% kh2=10; Nr2=0; reshuffle_flag2='none'; smooth_flag2='mean'; GRtype_flag2='dM'; NLEtype_flag2='comp';
% Generally B10dm reduces the confidence in bound assessments.
% Tends to improve on the failures of B10 with sporadic unbound cases.
% Is less certain about quantifiying between similar Mmax model types.
% 



