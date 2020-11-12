f_samp = 330e3;
f_trans=4e3;
%Band Edge specifications
fs1 = 75.9e3;
fp1 = 79.9e3;
fp2 = 99.9e3;
fs2 = 103.9e3;

%Defining normalized Cutoff frequencies 
Wc1 = (fp1-f_trans/2)*2*pi/f_samp
Wc2  = (fp2+f_trans/2)*2*pi/f_samp

%Defining Kaiser parameters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

dw= f_trans*2*pi/f_samp;                        %finding minimum transition width
N_min = 1+ ceil((A-8) / (2.285*dw)) ;         %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+17;                                 %final value found using trial and error

%Defining ideal BPF impulse response of length "n" and difference of 2 LPF
bp_ideal = ideal_lp(Wc2,n) - ideal_lp(Wc1,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
%which is obtained as rectangular in this case
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win               %multiplying the window and the ideal BPF response

%Plotting Magnitude (in dB) and Phase Responses in terms of normalized frequencies
fvtool(FIR_BandPass,'freq');    
figure(4)

%Plotting the time domain sequence
N=[-(n-1)/2:1:(n-1)/2]
stem(N,FIR_BandPass,'filled','Markersize',3)
xlim([-33,33])
ylim([-0.15,0.25])
xlabel('n')
ylabel('h[n]')
grid

hold off
%Magnitude response for un-normalized frequencies
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
figure(2)
plot(f/1000,abs(H));
%plot(-f/1000,abs(H),'b');            %for negative frequency axis response            
yline(0.15,'r');                    %marking important values of response and frequency
xline(79.9,'g--');
xline(75.9,'m--');
xline(99.9,'g--');
xline(103.9,'m--');
yline(0.85,'r');
yline(1.15,'r');
xlabel('Frequency (in kHz)')
ylabel('Magnitude Response')
legend('Magnitude Response','tolerances','Passband Edge','Stopband Edge')
grid