clc;clear;close all;
file_path='.\samples\Autel_evoII_pro_uplink_downlink_fs100e6.std';
fs=100e6;
fc=2440e6;


fid = fopen(file_path, 'r');
if fid == -1, error('无法打开文件: %s', file_path); end
fseek(fid, 50, 'bof');                          % pass 50 bytes 
rdata = fread(fid, inf, 'int16');
fclose(fid);
iq_data = rdata(1:2:end) + 1j * rdata(2:2:end);



%%
% downlink data

close all

downlinkdata = iq_data(160000:270000);

figure();plot(abs(downlinkdata))

pwelch(downlinkdata,[],[],[],fs,'centered');

%
offset = 16.4948e6;
downlinkdata=downlinkdata.*exp(1i*2*pi*(-offset)/fs*(0:length(downlinkdata)-1).');
pwelch(downlinkdata,[],[],[],fs,'centered');

 % stft
window = hamming(128); % window
noverlap = 0;        
nfft = 4096;           
figure,spectrogram(downlinkdata, window, noverlap, nfft, fs, 'yaxis','centered');


% xcorr
AA=xcorr(downlinkdata); 
figure,plot(abs(AA));



% 96668-90001 =  6667
% fs:  100MHz  6667点  bandwidth  9.2M
% 9.2/100*6666 
% 15.36MHz   1024     delta 15khz
% 15.36MHz/1024
% OFDM tu  = 1/Δf=1/15k≈66.67μs，

pa = fir1(64,1.2*15.36e6/fs);
downlinkdata= filter(pa,1,downlinkdata);
[num den] = rat(15.36e6/fs);
downlinkdata = resample(downlinkdata,num,den);  %resampling 
figure,plot(abs(downlinkdata));


% xcorr
AA1=xcorr(downlinkdata); %  1024
figure,plot(abs(AA1));


downlinkdata = downlinkdata./max(abs(downlinkdata));

cplen = 1024/16;
for i=1:length(downlinkdata)-1024-cplen
    C(i)=sum(downlinkdata(i:i+cplen).*conj(downlinkdata(i+1024:i+cplen+1024)));
end
figure,plot(abs(C));

the= 4;
value = abs(C(1,:)).*(abs(C(1,:))>the);
pks = [];
locs = [];
[pks(1,:), locs(1,:)] = findpeaks(value,'minpeakdistance', 1000);%1020

ind0=locs(1,1);
cplen = 1024/16-1;
figure,plot(angle(downlinkdata(ind0:ind0+cplen)));hold on,plot(angle(downlinkdata(ind0+1024:ind0+cplen+1024)));

df1=angle(sum(downlinkdata(ind0:ind0+cplen).*conj(downlinkdata(ind0+1024:ind0+cplen+1024))))*15.36e6/1024/2/pi;
rx_signal_fine=downlinkdata.*exp(1i*2*pi*(df1)/(15.36e6)*(0:length(downlinkdata)-1).');

figure,plot(angle(rx_signal_fine(ind0:ind0+cplen)));hold on,plot(angle(rx_signal_fine(ind0+1024:ind0+cplen+1024)));

rx_signal_fine = rx_signal_fine./max(abs(rx_signal_fine));


for loop = 1: length(locs)
    ofdm(:,loop) = rx_signal_fine(locs(loop):locs(loop)+1024+1024/16-1);
end


 sym = fftshift(fft(ofdm(1+1024/16:1024+1024/16,:),1024));

 figure();plot(abs(sym(:,1)))
 
spe = fft(sym(:,2)).^2,1024*8;

 figure();plot(abs(spe))
 
 for lll = 1:5
    scatterplot(sym(512-8:512+8,lll)./max(abs(sym(:,lll))));
 end
 
 
 
 
 
 
 