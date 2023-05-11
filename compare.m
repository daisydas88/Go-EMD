clc;
close all;
clear;

data = importdata('preictal50.mat');
data = sgolayfilt(data', 5, 33);
zz = data;
imf = eemd(data, 12, 1, 1.2);
[nrow, ncol] = size(imf);

featurematrix = [];

%code for PSD for table II
Fs=200;
for ii=1:nrow
    
    p=imf(ii, :)';
    [Pxx, F] = periodogram(p, rectwin(length(p)), length(p), Fs);
%     delta = bandpower(Pxx,F,[0 4],'psd');
%     theta = bandpower(Pxx,F,[4 8],'psd');
    alpha = bandpower(Pxx,F,[8 12],'psd');
    beta = bandpower(Pxx,F,[12 35],'psd');
    gamma = bandpower(Pxx,F,[35 70],'psd'); 

%     featurematrix(ii, 1) = delta;
%     featurematrix(ii, 2) = theta;
    featurematrix(ii, 1) = alpha;
    featurematrix(ii, 2) = beta;
    featurematrix(ii, 3) = gamma;


end


%% No. of iteration = goal * ens * itr [from emd.m]
%% 

rs = [];
for ii=1:ncol
    rs(1, ii) = sum(imf(:, ii));
end
c=data;
figure
subplot(2, 1, 1);
plot(data);
subplot(2, 1, 2);
plot(rs);
%code for variance of whole signal
wholevar = var(c)
%code for variance of imfs of signal after decomposition table IV
m=imf(:,1)';
y1 = var(m)
m=imf(:,2)';
y2 = var(m)
m=imf(:,3)';
y3 = var(m)
m=imf(:,4)';
y4 = var(m)
m=imf(:,5)';
y5 = var(m)
m=imf(:,6)';
y6 = var(m)
m=imf(:,7)';
y7 = var(m)
m=imf(:,8)';
y8 = var(m)
m=imf(:,9)';
y9 = var(m)
m=imf(:,10)';
y10 = var(m)
m=imf(:,11)';
y11 = var(m)
m=imf(:,12)';
y12 = var(m)
ivar=y1+y2+y3+y4+y5+y6+y7+y8+y9+y10+y11+y12;
%mean
wholemean = mean(c)
%code for mean of imfs of signal after decomposition for table IV
n=imf(:,1)';
my1 = mean(n)
n=imf(:,2)';
my2 = mean(n)
n=imf(:,3)';
my3 = mean(n)
n=imf(:,4)';
my4 = mean(n)
n=imf(:,5)';
my5 = mean(n)
n=imf(:,6)';
my6 = mean(n)
n=imf(:,7)';
my7 = mean(n)
n=imf(:,8)';
my8 = mean(n)
n=imf(:,9)';
my9 = mean(n)
n=imf(:,10)';
my10 = mean(n)
n=imf(:,11)';
my11 = mean(n)
n=imf(:,12)';
my12 = mean(n)
imean=my1+my2+my3+my4+my5+my6+my7+my8+my9+my10+my11+my12;
%code for standard deviation of whole signal
wholestd=std(c)
%code for standard deviation of imfs for tableIV
ss=imf(:,1)';
sv1=std(ss)

ss=imf(:,2)';
sv2=std(ss)
ss=imf(:,3)';
sv3=std(ss)
ss=imf(:,4)';
sv4=std(ss)
ss=imf(:,5)';
sv5=std(ss)
ss=imf(:,6)';
sv6=std(ss)
ss=imf(:,7)';
sv7=std(ss)
ss=imf(:,8)';
sv8=std(ss)
ss=imf(:,9)';
sv9=std(ss)
ss=imf(:,10)';
sv10=std(ss)
ss=imf(:,11)';
sv11=std(ss)
ss=imf(:,12)';
sv12=std(ss)
isd=sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10+sv11+sv12;
%code for skewness of whole signal and resconstructed signal for table IV
wholesk = skewness(c)
ss=imf(:,1)';
sk1=skewness(ss)

ss=imf(:,2)';
sk2=skewness(ss)
ss=imf(:,3)';
sk3=skewness(ss)
ss=imf(:,4)';
sk4=skewness(ss)
ss=imf(:,5)';
sk5=skewness(ss)
ss=imf(:,6)';
sk6=skewness(ss)
ss=imf(:,7)';
sk7=skewness(ss)
ss=imf(:,8)';
sk8=skewness(ss)
ss=imf(:,9)';
sk9=skewness(ss)
ss=imf(:,10)';
sk10=skewness(ss)
ss=imf(:,11)';
sk11=skewness(ss)
ss=imf(:,12)';
sk12=skewness(ss)

isk=sk1+sk2+sk3+sk4+sk5+sk6+sk7+sk8+sk9+sk10+sk11+sk12;
%CODE for IMF plot Figure 10
 for ii=1:12
  subplot(12,1,ii);
  plot(imf(ii,:));
  xlabel('Time');
ylabel('Fq');
  title(sprintf("IMF: %d",ii));
if ii==12
   subplot(12,1,ii);
  plot(imf(ii,:));
  xlabel('Time');
ylabel('Fq');
  title(sprintf("RES")); end
 end
%code for reconstruction figure 8
%   rs=imf(:,1)'+imf(:,2)'+imf(:,3)'+imf(:,4)'+imf(:,5)'+imf(:,6)'+imf(:,7)'+imf(:,8)'+imf(:,9)'+imf(:,10)'+imf(:,11)'+imf(:,12)'
rs=imf(1,:)+imf(2,:)+imf(3,:)+imf(4,:)+imf(5,:)+imf(6,:)+imf(7,:)+imf(8,:)+imf(9,:)+imf(10,:)+imf(11,:)+imf(12,:)
yy=data;
figure;
subplot(2,1,1);
plot(rs(1:999));
xlabel('time(s)')
title('Reconstructed signal from EEMD')
subplot(2,1,2);
plot(yy(1:999));
xlabel('time(s)')
ylabel('freq')
title('Original Signal')

%code for SNR
 r = snr(zz, rs)

 %code for RMSE
 E = rmse(rs,data)


 %nmse
 n = calNMSE(zz, rs);


 % SSIM
ss = ssim(zz, rs);
 % CC
cc = corrcoef(zz, rs);

for ii=1:12
%Code for figure 6
  figure;
  img = plot(imf(ii, 1:1000));
  xlabel('Time');
  ylabel('Freq');
  title(strcat('IMF ', int2str(ii)));
  outfile = strcat('EEMD_IMF_', int2str(ii), '.png')
  saveas(img, outfile);
end


  figure;
  fig = plot(imf(13,1:1000));
  xlabel('Time');
  ylabel('Freq');
  title('RES');
  saveas(fig, 'EEMD_RES.png')
  %code for reconstruction signal in one plot

yyy=zz';
figure;
plot(rs(1:999),'-b'); hold on;
plot(yyy(1:999),'-r');
ylabel('freq')
xlabel('time(s)')
title('Reconstructed signal from EEMD') ;
legend('Reconstructed Signal by EEMD','Original EEG signal')
hold off;