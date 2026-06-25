clc;clear;close all;
file_path='.\samples\Autel_evoII_pro_uplink_downlink_fs100e6.std';
fs=100e6;
fc=2440e6;


fid = fopen(file_path, 'r');
if fid == -1, error('无法打开文件: %s', file_path); end
fseek(fid, 50, 'bof');                          % 跳过50字节STD头
rdata = fread(fid, inf, 'int16');
fclose(fid);
iq_data = rdata(1:2:end) + 1j * rdata(2:2:end);



%%
% uplink data
uplinkdata = iq_data(70000:160000);

pwelch(uplinkdata,[],[],[],fs,'centered');

%
offset = -21.4966e6;
uplinkdata=uplinkdata.*exp(1i*2*pi*(-offset)/fs*(0:length(uplinkdata)-1).');
pwelch(uplinkdata,[],[],[],fs,'centered');

 % 时频图
window = hamming(128); % 窗函数
noverlap = 0;        % 重叠样本数
nfft = 4096;           % FFT点数
figure,spectrogram(uplinkdata, window, noverlap, nfft, fs, 'yaxis','centered');


% 自相关有效符号估计 
AA=xcorr(uplinkdata); 
figure,plot(abs(AA));
% 96672-90001 =  6671

% bw = 1Mhz

 
 
 
 
 