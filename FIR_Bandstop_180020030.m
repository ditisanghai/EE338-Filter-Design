f_samp = 260e3;
f_trans=4e3;
%Band Edge specifications
fs1 = 61.9e3;
fp1 = 57.9e3;
fp2 = 85.9e3;
fs2 = 81.9e3;

%Finding Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

dw= 4*2*pi/260  ;                        %minimum transition bandwidth
N_min = ceil((A-8) / (2.285*dw));     %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+16                          %Obtained using trial and error

%Ideal BSF impulse response of length using multiple LPF

bs_ideal =  ideal_lp(pi,n) -ideal_lp((fp2-2000)*pi*2/f_samp,n) + ideal_lp((fp1+2000)*pi*2/f_samp,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
%which turns out to be a rectangular window in this case
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win           %multiplying the window and the ideal BPF response

%Plotting Magnitude (in dB) and Phase Responses in terms of normalized frequencies
fvtool(FIR_BandStop,'freq'); 

%Plotting the time domain sequence
figure(4)
N=[-(n-1)/2:1:(n-1)/2]
stem(N,FIR_BandStop,'filled','Markersize',3)
xlim([-27,27])
ylim([-0.15,0.85])
xlabel('n')
ylabel('h[n]')
grid


%Magnitude response for un-normalized frequencies
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
figure(3)
plot(f/1000,abs(H));
%plot(-f/1000,abs(H),'b');            %for negative frequency axis response
yline(0.85,'r');                        %marking important values of response and frequency
xline(57.9,'g--');
xline(61.9,'m--');
xline(81.9,'m--');
xline(85.9,'g--');
yline(0.15,'r');
yline(1.15,'r');
xlabel('Frequency (in kHz)')
ylabel('Magnitude Response')
legend('Magnitude Response','tolerances','Passband Edge','Stopband Edge','location','southwest')
grid

