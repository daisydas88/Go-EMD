% For reconstructed spectral features
% code for peak frequency
fs=200;
k= fft(s);

f=(0:length(k)-1)*fs/length(k);
plot(f,abs(k));
pow= k.*conj(k);
meanFreq = meanfreq(k);
medianFreq= medfreq(k);
totalPower = sum(pow);
[maxAmp, maxAmpIdx] = max(abs(k));
peakFreq = f(maxAmpIdx);
disp(peakFreq);
% code for no of peaks
% rs is the reconstructed signal
a= findpeaks(rs);
disp(length(a));
% code for peak amplitude
pa=max(abs(rs));
disp(pa);