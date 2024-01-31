%% 设置参数
clear;clc;
Nk = 32; % 子载波个数
Nfft = 32; % fft长度
Nframe = 1; % 一帧中有几个OFDM符号
M = 4; % 调制符号所含比特（改为4：4QAM，16：16QAM,32：32:QAM）
SR = 250000; % 符号速率
BR = SR * log2(M); % 比特率（根据调制方式计算比特率）
NGI = 32; % 保护间隔长度
EbN0s = 0:2:10; % 信噪比
Nsym = Nfft + NGI; % 系统长度
bers = zeros(1, length(EbN0s));

PowerTDL_dB = [0 -8 -17 -21 -25];   % TDL中信道抽头的功率,dB为单位
Delay = [0 3 5 6 8];                % TDL中信道时延
PowerTDL = 10.^(PowerTDL_dB/10);    % TDL中信道抽头的功率
Nchannel=length(PowerTDL_dB);       % 信道抽头数
Tau_maxTDL = Delay(end)+1;          % 最大时延除以帧长,就是归一化的最大时延
fprintf('EbN0 \t \t ber\t\t\t per\t\t\t nloop \t\t \n');
%% 函数主体

for kk = 1:length(EbN0s)
    % rng('default')          % 初始化随机种子
    EbN0 = EbN0s(kk);
    nloop = 10000; % 发送多少帧
    n_biterror = 0; % 错误的数据
    n_bitdata = 0; % 一共发送了多少数据
    n_packeterror = 0; % 有多少错误帧
    n_packetdata = 0; % 发送了多少帧

    for ii = 1:nloop
        % 生成一帧数据，串并转换，并QPSK，生成一帧
        frame_FDserial = randi([0 1], 1, Nk * Nframe * log2(M));% 发送的是bit（根据M修改生成的比特数）
        frame_FDparallel = reshape(frame_FDserial, Nk, Nframe * log2(M)); % 串并转换
        frame_mod = QAMMod(frame_FDparallel, Nk, Nframe, M); %调制（改为QAM调制）
        % IFFT
        power_FT = sum(abs(frame_mod(:)).^2) / (Nk * Nframe);% 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_mod); % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft); % ifft
        power_TD = sum(sum(abs(frame_ifft).^2)) / Nk / Nframe; % 计算下IFFT前的能量，DT表示时域
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP"); % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI, 1, Nsym * Nframe* log2(M));
        % Channel         
        
        channel = Rayleigh_model(Nchannel, PowerTDL);
        h = zeros(1, Tau_maxTDL);
        h(Delay+1) = channel;
        frame_conv = conv(frame_TDserial, h);
        frame_fading = frame_conv(:,1:length(frame_TDserial));        % 看似是线性卷积，实际上由于CP变成了循环卷积
        % 添加高斯白噪声
        power_TDserial = sum(abs(frame_TDserial).^2)/Nk/Nframe;     % 计算出的能量和理论不符啊，没发现问题在哪
        EsN0 = EbN0 + 10*log10(M);                                  % 根据信噪比计算噪声能量，幅值，然后加在信号上
        N0 = power_TDserial .* 10.^(-EsN0/10);
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
        frame_recieved = frame_fading + noise_msg;
        
        
        
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved, Nsym, Nframe* log2(M));
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        
        % 信道均衡
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft));
        frame_equalization = frame_recieved_FD ./ repmat(H, 1, Nframe* log2(M));
        % QPSK解调
        frame_demod = QAMDemod(frame_equalization, Nk, Nframe, M); %改为QAM解调
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk * Nframe * log2(M)); %修改输出比特数

        % 计算error
        n_biterror_tmp = sum(abs(frame_output - frame_FDserial));
        n_bitdata_tmp = length(frame_FDserial);
        n_biterror = n_biterror + n_biterror_tmp;
        n_bitdata = n_bitdata + n_bitdata_tmp;

        if n_biterror_tmp ~= 0
            n_packeterror = n_packeterror + 1;
        end

        n_packetdata = n_packetdata + 1;
    end

    % 计算在当前信噪比下的误码率
    per = n_packeterror / n_packetdata;
    ber = n_biterror / n_bitdata;
    bers(kk) = ber;
    fprintf('%f\t%e\t%e\t%d\t\n', EbN0, ber, per, nloop);
end

semilogy(EbN0s, bers, '-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend(strcat(num2str(M),'QAM实验曲线'));

function outs = QAMMod(input_data,nk,nframe,M)
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*log2(M))
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% M: 调制数
% 输出
% outs: (nk,nframe),输出a+bi
% out_coss:(nk,nframe)，输出实部
% out_sins:(nk,nframe)，输出虚部

if nargin < 4                   % 设置默认值
    M = 4; % 默认为4QAM
end

% 将输入二进制数据映射成QAM符号
symbols = qammod(input_data, M);

% 将QAM符号按照OFDM格式进行串并转换
outs = reshape(symbols, nk, nframe* log2(M));
end

function outputs = QAMDemod(input_data,nk,nframe,M)
% 输入
% input_data: (Nk, Nframe), 一个频域的复数，会被拆开解调
% nk: 频域并联
% nframe: 一帧包含符号数
% M: 调制数
% 输出
% outputs：(Nk, Nframe*log2(M)), 解调后的比特流

if nargin < 4                   % 设置默认值
    M = 4; % 默认为4QAM
end

% 将输入QAM符号按照OFDM格式进行并串转换
symbols = input_data(:);

% 将QAM符号解调为二进制比特流
outputs = qamdemod(symbols, M);
end

function output_TD = AddGI(input_TD, nfft, nGI, nframe, type_GI)
if type_GI=="CP"    % 实现CP
    output_TD = [input_TD(nfft-nGI+1:nfft, :); input_TD(1:nfft, :)];
elseif type_GI=="ZP" % 实现ZP
    output_TD = [zeros(nGI,nframe); input_TD(1:nfft, :)];
end
end
function output_TD = RemoveGI(input_TD,nfft,nGI)
% 输入
% input_TD: (Nsym,Nframe)输入的并联时域数据
% nfft：fft长度
% nGI: GI长度
% 输出
% output_TD: (Nfft,Nframe)去掉GI后的并联时域数据
    output_TD = input_TD(nGI+1:nfft+nGI,:);
end
function H=Rayleigh_model(nchannel, power_channel)
% 瑞利衰落信道
% 输入
% nchannel： 多径信道的个数
% power_channel：（1, nchannel），每一个信道的功率
% 输出
% H:(1, nchannel),一个瑞利信道,符合高斯分布的nchannel个随机数，代表着衰落
H = (randn(1,nchannel)+1j*randn(1,nchannel)).*sqrt(power_channel/2);
% 功率除以二的原因是瑞利分布的E(x^2)=2\sigma^2
end