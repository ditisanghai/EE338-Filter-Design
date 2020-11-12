%Defining Butterworth Analog LPF parameters
Wc = 1.071;              %Cut-off frequency
N = 7;                   %Minimum order 

%Calculating poles of Butterworth polynomial of degree N=7 
p1 = Wc*cos(pi/2 + pi/14) + i*Wc*sin(pi/2 + pi/14);
p2 = Wc*cos(pi/2 + pi/14) - i*Wc*sin(pi/2 + pi/14);
p3 = Wc*cos(pi/2 + pi/14+pi/7) + i*Wc*sin(pi/2 + pi/14+pi/7);
p4 = Wc*cos(pi/2 + pi/14+pi/7) - i*Wc*sin(pi/2 + pi/14+pi/7);
p5 = Wc*cos(pi/2 + pi/14+2*pi/7) + i*Wc*sin(pi/2 + pi/14+2*pi/7);
p6 = Wc*cos(pi/2 + pi/14+2*pi/7) - i*Wc*sin(pi/2 + pi/14+2*pi/7);
p7 = Wc*cos(pi/2 + pi/14+3*pi/7) + i*Wc*sin(pi/2 + pi/14+3*pi/7);

%Band Edge speifications in kHz
fp1 = 79.9;
fs1 = 75.9;
fs2 = 103.9;
fp2 = 99.9;

%Transforming Band Edge specifications using Bilinear Transformation
f_samp = 330;       %in kHz  
wp1 = tan(fp1*pi/f_samp);
ws1 = tan(fs1*pi/f_samp);
ws2 = tan(fs2*pi/f_samp);
wp2 = tan(fp2*pi/f_samp);

%Defining Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7],Wc^N);   %transfer function with poles p1-p7 and numerator Wc^N
                                                 

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);       %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));       %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation in z
discrete_lpf(z)= analog_lpf((z-1)/(z+1));



%Finding coefficients of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to get coefficients
ns = sym2poly(expand(ns))  ;                       
ds = sym2poly(expand(ds))                           %getting coefficients in decimal form
k = ds(1);                                         % setting coefficients of highest degree of denominator to 1
ds = ds/k;                                          %final coefficients of the numerator and denominator polynomial
ns = ns/k;

%Finding coefficients of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to get coefficients                   
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %getting coefficients in decimal form
k = dz(1);                                             %setting coefficients of highest degree of denominator to 1
dz = dz/k;                                              %final coefficients of the numerator and denominator polynomial
nz = nz/k;

%Plotting Magnitude (in dB) and Phase Responses in terms of normalized frequencies 
fvtool(nz,dz,'freq');                                          


%Magnitude response for un-normalized frequencies
[H,f] = freqz(nz,dz,1024*1024, f_samp*1e3);
plot(f/1000,abs(H),'b');
%hold on
%plot(-f/1000,abs(H),'b');            %for negative frequency axis response                
grid
xlabel('Frequency (in kHz)');
ylabel('Magnitude Response');
xline(79.9,'--m');                      %marking important values of response and frequency on the plot
xline(75.9,'--g');
yline(0.15,'r');
xline(99.9,'--m');
xline(103.9,'--g');
yline(1.15,'r');
yline(0.85,'r');
legend('Magnitude Response','Passband edge','Stopband edge','Tolerances','location','northwest');
