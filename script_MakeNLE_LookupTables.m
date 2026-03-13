% Script to make log-likelihood lookup tables for dMlrg(i).
% This is my workaround for the slow nLL computation, since the dM convolution integral prevents code vectorization.
clear;

% Define some constants.
Nc=1e4;
m1=0.0;
dx=1e-4;
d1=0.01;
d2=0.1;
d3=0.01;
kd=10;
type='norm';
O=-1;

% Define the relevant axes.
mx=0:dx:10;
ma=0:d1:5;
ba=d2:d2:3.00;
m2a=d3:d3:5;


% Derive some more variables.
a=log10(Nc);

% Preallocate.
Sa=struct('D',zeros([length(ma) length(ba) length(m2a)]),'GI',[]);

% Loop over all the order statistics, b-values, and Mmax.
for i=1:kd
    i
    for j=1:length(ba)
        for l=1:length(m2a)
            
            % Compute the PDF for these values, as a fxn of the magnitude axis.
            [dMa,dMLpdf,dMLcdf,dMLsvf]=GR_NLE(mx, m1,m2a(l), a,ba(j), i,type,O);
            PDF=interp1(dMa,dMLpdf,ma,'linear')';
            PDF(ma>=m2a(l))=1e-16;
            
            % Stuff values into the structure/table.
            Sa(i).D(:,j,l)=PDF;
        end
    end
end

% Make the gridded interpolant.
[x,y,z]=ndgrid(ma,ba,m2a);
for i=1:length(Sa)
    Sa(i).GI=(griddedInterpolant(x,y,z,log(Sa(i).D), 'linear'));
end

% Save the tables.
save('NLEtable_temp.mat','Sa','ma','ba','m2a');

% Make a smaller version of the table and save that too.
for i=1:length(Sa)
    Sa(i).D=[];
end
save('NLEtable_small_temp.mat','Sa','ma','ba','m2a');




%%
% Compare the table interpolation against the numerical function.
it=2;
bt=1.00;
m2t=10;
mt=0:0.005:min([max(m2t),10]);
load('NLEtable.mat','Sa','ma','ba','m2a');

% Get the PDF from the lookup table and the actual function.
%PDFt=(interp3(ba,ma,m2a,(S(it).D), bt,mt,m2t, 'linear'));
PDFt=exp(Sa(it).GI(mt,bt*ones(size(mt)),5*ones(size(mt))));
[dM,dMLpdf,dMLcdf,dMLsvf]=GR_NLE(0:1e-4:10, 0,m2t(1), 0,bt(1), it,'norm',-1);

% Plot.
figure(1); clf;
subplot(211);
semilogy(dM,dMLpdf,'-b'); hold on;
semilogy(mt,(PDFt),'xr');
xlabel('\DeltaM'); ylabel('PDF');
subplot(212);
y=interp1(dM,dMLpdf,mt,'linear');
plot(mt,100*((PDFt)-y)./y,'-b');
xlabel('\DeltaM'); ylabel('PDF relative residual (%)');

% Compare LL values.
sum(log(PDFt))
GR_NLE_LLtable(mt,m1,m2t*ones(size(mt)),a,bt*ones(size(mt)),it, Sa,ma,ba,m2a)


