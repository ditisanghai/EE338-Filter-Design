%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;                      %Minimum order of the filter

%Poles of the Chebyshev Polynomial of order 4
p1=-0.294921 + 0.401709i;
p2=-0.294921-0.401709i;
p3=-0.12216 + 0.96981i;
p4=-0.12216-0.96981i;        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))] ;       % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Defining Band Edge speifications in kHz
fs1 = 61.9;
fp1 = 57.9;
fp2 = 85.9;
fs2 = 81.9;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;                       %in kHz
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog LPF transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s +W0*W0));     %bandpass transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation in z

%Finding coefficients of analog BSF
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to get coefficients
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %getting coefficients in decimal form
k = ds(1);                                          %setting coefficients of highest degree of denominator to 1
ds = ds/k ;                                          %final coefficients of the numerator and denominator polynomial
ns = ns/k;

%Finding the coefficients of discrete BSF
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification to get coefficients
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %getting coefficients in decimal form
k = dz(1);                                          %final coefficients of the numerator and denominator polynomial
dz = dz/k;
nz = nz/k;
%Plotting Magnitude (in dB) and Phase Responses in terms of normalized frequencies
fvtool(nz,dz,'freq');                                      %frequency response in dB

%Magnitude response for un-normalized frequencies
[H,f] = freqz(nz,dz,1024*1024, f_samp*1e3);
plot(f/1000,abs(H));                                    
%plot(-f/1000,abs(H),'b');            %for negative frequency axis response
xline(61.9,'--m');                      %marking important values of response and frequency on the plot
xline(57.9,'--g');
yline(0.15,'r');
xline(81.9,'--m');
xline(85.9,'--g');
yline(1.15,'r');
yline(0.85,'r');
xlim([0,130]);
ylim([0,1.25]);
legend('Magnitude Response','Stopband edge','Passband edge','Tolerances','location','southwest');
xlabel('Frequency (kHz)');
ylabel('Magnitude Response');
grid