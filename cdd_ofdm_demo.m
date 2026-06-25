%% CDD-OFDM 2Tx-1Rx 仿真 MATLAB R2021b+ 无维度报错
clear; clc; close all;
%% ====================== 系统参数 ======================
N_fft   = 1024;      % FFT点数
N_cp    = 64;       % CP长度
N_sym   = 200;      % 符号数(减少加快仿真)
mod_order = 2;      % QPSK
M_tx    = 2;        % 2发射天线
delta   = [0,8];    % 循环延迟点数
EbN0_dB = 30;       % 信噪比

%% 梳状导频
pilot_interval = 8;
pilot_pos = (1:pilot_interval:N_fft).';  % 列向量统一维度
N_pilot = length(pilot_pos);
pilot_seq = (1-2*randi([0,1],N_pilot,1)) + 1j*(1-2*randi([0,1],N_pilot,1));

%% ====================== 发送端 CDD-OFDM ======================
bit_total = N_sym * (N_fft - N_pilot) * mod_order;
tx_bits = randi([0,1], bit_total, 1);
bit_idx = 1;

rx_time_wave = zeros(N_sym*(N_fft+N_cp), 1);

for sym_idx = 1:N_sym
    X_freq = zeros(N_fft, 1);  % 统一为列向量 Nfft×1
    data_pos = setdiff((1:N_fft).', pilot_pos);
    N_data = length(data_pos);
    
    bit_pack = tx_bits(bit_idx : bit_idx+N_data*mod_order-1);
    bit_idx = bit_idx + N_data*mod_order;
    
    % QPSK映射
    bit_mat = reshape(bit_pack, mod_order, []).';
    sym_I = 1 - 2*bit_mat(:,1);
    sym_Q = 1 - 2*bit_mat(:,2);
    X_data = (sym_I + 1j*sym_Q)/sqrt(2);
    
    X_freq(data_pos) = X_data;
    X_freq(pilot_pos) = pilot_seq;
    
    % IFFT 基础时域符号
    x_base = ifft(X_freq, N_fft) * sqrt(N_fft);
    
    % CDD 多天线循环移位
    tx_ant = zeros(M_tx, N_fft);
    for ant = 1:M_tx
        d = delta(ant);
        tx_ant(ant,:) = circshift(x_base.', d);
    end
    
    % 加CP
    tx_ant_cp = zeros(M_tx, N_fft+N_cp);
    tx_ant_cp(:,1:N_cp) = tx_ant(:,end-N_cp+1:end);
    tx_ant_cp(:,N_cp+1:end) = tx_ant;
    
    % 瑞利多径信道 comm.RayleighChannel
    fs = 1e6;
    pathDelays = [0 1e-6 2.5e-6];
    avgPathGains = [0 -3 -6];
    
    rayChan1 = comm.RayleighChannel(...
        'SampleRate',fs,...
        'PathDelays',pathDelays,...
        'AveragePathGains',avgPathGains,...
        'MaximumDopplerShift',1);
    rayChan2 = comm.RayleighChannel(...
        'SampleRate',fs,...
        'PathDelays',pathDelays,...
        'AveragePathGains',avgPathGains,...
        'MaximumDopplerShift',1);
    
    r1 = rayChan1(tx_ant_cp(1,:).');
    r2 = rayChan2(tx_ant_cp(2,:).');
    r_sym = r1 + r2;
    
    rx_time_wave((sym_idx-1)*(N_fft+N_cp)+1 : sym_idx*(N_fft+N_cp)) = r_sym;
end

%% ====================== 加复高斯噪声 ======================
Es = mean(abs(rx_time_wave).^2);
Eb = Es / mod_order;
sigma_n = sqrt(Eb / (10^(EbN0_dB/10)));
rx_noisy = rx_time_wave + sigma_n * (randn(size(rx_time_wave)) + 1j*randn(size(rx_time_wave)))/sqrt(2);
noise_var = sigma_n^2;

pwelch(rx_time_wave,[],[],[],fs,'centered');


%% ====================== CDD 解调核心（修复维度） ======================
rx_out_bits = [];
H_eq_last = zeros(N_fft,1);

for sym_idx = 1:N_sym
    sym_start = (sym_idx-1)*(N_fft+N_cp) + 1;
    r_full = rx_noisy(sym_start : sym_start + N_fft + N_cp - 1);
    r_nocp = r_full(N_cp+1 : end);
    
    % FFT 输出 Nfft×1 列向量
    R_freq = fft(r_nocp, N_fft) / sqrt(N_fft);
    
    % LS信道估计
    R_pilot = R_freq(pilot_pos);
    H_pilot_ls = R_pilot ./ pilot_seq;
    
    % 样条插值输出 H_eq 列向量 Nfft×1
    H_eq = spline(pilot_pos, H_pilot_ls, (1:N_fft).');
    H_eq_last = H_eq;
    
    % ========== 关键修复：全部同维度(Nfft×1)，无维度冲突 ==========
    H_conj = conj(H_eq);
    denom = abs(H_eq).^2 + noise_var;
    X_hat = (H_conj ./ denom) .* R_freq;
    
    % 提取数据子载波
    X_data_hat = X_hat(data_pos);
    
    % QPSK硬解调
    bit_I = real(X_data_hat) < 0;
    bit_Q = imag(X_data_hat) < 0;
    demod_bits = [bit_I, bit_Q].';
    demod_bits = demod_bits(:);
    rx_out_bits = [rx_out_bits; demod_bits];
end

%% 误码统计
bit_err = sum(xor(tx_bits, rx_out_bits));
BER = bit_err / bit_total;
fprintf('EbN0 = %d dB, BER = %.4e\n', EbN0_dB, BER);

%% 绘制CDD等效信道
figure('Color','w');
plot(1:N_fft, 20*log10(abs(H_eq_last)));
title('CDD等效频域信道幅度 |H_{eq}[k]|');
xlabel('子载波索引 k');
ylabel('幅度(dB)');
grid on;