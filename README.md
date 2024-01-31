# 基于MATLAB实现的OFDM仿真调制解调，BPSK、QPSK、4QAM、16QAM、32QAM，加性高斯白噪声信道、TDL瑞利衰落信道

# 相关链接

[OFDM中的帧（frame）、符号（symbol）、子载波（subcarriers）、导频（Pilot）、保护间隔（guard）的关系图解以及代码详解--MATLAB](http://t.csdnimg.cn/OIeEU)

[【学习笔记】OFDM的原理和技术介绍以及仿真结果分析附代码--MATLAB](http://t.csdnimg.cn/78aAz)



# 1 加性高斯白噪声信道

在本程序中，通过MATLAB仿真了OFDM的发射、信道、接收解调的过程。支持的BPSK、QPSK、多种QAM的解调方式，并计算了不同信噪比下的误比特率。高斯白噪声信道下，不需要信道估计和信道均衡，但是在衰落信道中必须要信道估计或信道均衡，或者两个都使用，才能正确的解调。以下在衰落信道中的解调只实现了信道均衡。

## 1.1 BPSK和QSPK

```matlab
%-----------------------仿真OFDM---------------------------%
%% 设置参数
clear;clc;
Nk = 128; % 子载波个数
Nfft = 128; % fft长度
Nframe = 6; % 一帧中有几个OFDM符号

M = 2; % 调制符号所含比特（改为1：BPSK，2：QPSK）
SR = 250000; % 符号速率
BR = SR .* M; % 比特率
NGI = 32; % 保护间隔长度
EbN0s = 0:1:12; % 信噪比
Nsym = Nfft + NGI; % 系统长度
bers = zeros(1, length(EbN0s));
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
        frame_FDserial = randi([0 1], 1, Nk * Nframe * M);% 发送的是bit
        frame_FDparallel = reshape(frame_FDserial, Nk, Nframe * M); % 串并转换
        if M==1
            frame_mod = BPSKMod(frame_FDparallel, Nk, Nframe); % BPSK调制
        elseif M==2
            frame_mod = QPSKMod(frame_FDparallel, Nk, Nframe); % QPSK调制
        end
        % IFFT
        power_FT = sum(abs(frame_mod(:)).^2) / (Nk * Nframe);% 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_mod); % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft); % ifft
        power_TD = sum(sum(abs(frame_ifft).^2)) / Nk / Nframe; % 计算下IFFT前的能量，DT表示时域
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP"); % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI, 1, Nsym * Nframe);
        x = 1:1:160;
        % 加性白高斯噪声信道         
        power_TDserial = sum(abs(frame_TDserial).^2) / length(frame_TDserial); % 计算发送序列的能量
        EsN0 = EbN0 + 10 * log10(M); % 根据信噪比计算Es/N0(dB)噪声能量，幅值，然后加在信号上
        N0 = power_TDserial / (10^(EsN0 / 10));
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1i * randn(size(frame_TDserial)));
        frame_recieved = frame_TDserial + noise_msg;
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved, Nsym, Nframe);
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        if M==1
            % BPSK解调
            frame_demod = BPSKDemod(frame_recieved_FD, Nk, Nframe);
        elseif M==2
            % QPSK解调
            frame_demod = QPSKDemod(frame_recieved_FD, Nk, Nframe);
        end
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk * Nframe * M);
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
if M==1
    legend('BPSK实验曲线');
else
    legend('QPSK实验曲线');
end

function outs = QPSKMod(input_data,nk,nframe,m)
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*M)
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% m: 调制数
% 输出
% outs: (nk,nframe),输出a+bi
% out_coss:(nk,nframe)，输出实部
% out_sins:(nk,nframe)，输出虚部

if nargin < 4                   % 设置默认值
    m = 2;
end

outs = zeros(nk,nframe);    % 初始化
out_cos = zeros(nk,nframe);
out_sin = zeros(nk,nframe);

input_data_pn = 2*input_data - 1;     % 把0，1映射到-1，1
A = 1/sqrt(2);                        % 归一化幅值

out_cos((1:nk),(1:nframe)) = input_data_pn((1:nk), (1:2:2*nframe-1)) .* A;    % 每次使用两列
out_sin((1:nk),(1:nframe)) = input_data_pn((1:nk), ((2:2:2*nframe))) .* A;

outs = out_cos + out_sin * 1j;
end

function outputs = QPSKDemod(input_data,nk,nframe)
% 输入
% input_data: (Nk, Nframe), 一个频域的复数，会被拆开解调
% nk: 频域并联
% nframe: 一帧包含符号数
% 输出
% outputs：(Nk, 2*Nframe), 解调后，多出一倍，全是01
outputs = zeros(nk,2*nframe);
A = 1/sqrt(2); 
input_data = input_data ./ A;
outputs((1:nk),(1:2:2*nframe-1)) = real(input_data((1:nk),(1:nframe)))>=0;
outputs((1:nk),(2:2:2*nframe)) = imag(input_data((1:nk),(1:nframe)))>=0;

end

function outs = BPSKMod(input_data,nk,nframe,m)
% BPSK调制函数
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*M)
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% m: 调制数，默认值为1
% 输出
% outs: (nk,nframe),输出复数
A = 1/sqrt(2);    % 归一化幅值
outs = input_data * 2 - 1;  % 将二进制0/1映射为复数-1/1
outs = A * outs;  % 乘以归一化幅值
end

function outputs = BPSKDemod(input_data,nk,nframe)
% BPSK解调函数
% 输入
% input_data: (Nk, Nframe), 一个频域的复数
% nk: 频域并联
% nframe: 一帧包含符号数
% 输出
% outputs：(Nk, Nframe), 解调后的二进制数据
outputs = real(input_data) >= 0;  % 判断实部是否大于等于0，大于等于0表示1，小于0表示0
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
```
![在这里插入图片描述](https://img-blog.csdnimg.cn/direct/b8f14041c34f463d8af29f8c22f275f1.png#pic_center)



## 1.2 QAM

M = 16; % 调制符号所含比特（改为4：4QAM，16：16QAM,32:32QAM）

```matlab
%% 设置参数
clear;clc;
Nk = 128; % 子载波个数
Nfft = 128; % fft长度
Nframe = 6; % 一帧中有几个OFDM符号
M = 16; % 调制符号所含比特（改为4：4QAM，16：16QAM,32:32QAM）
SR = 250000; % 符号速率
BR = SR * log2(M); % 比特率（根据调制方式计算比特率）
NGI = 32; % 保护间隔长度
EbN0s = 0:1:12; % 信噪比
Nsym = Nfft + NGI; % 系统长度
bers = zeros(1, length(EbN0s));
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
        x = 1:1:160;
        % Channel         
        power_TDserial = sum(abs(frame_TDserial).^2) / length(frame_TDserial); % 计算发送序列的能量
        EsN0 = EbN0 + 10 * log10(log2(M)); % 根据信噪比计算Es/N0(dB)噪声能量，幅值，然后加在信号上
        N0 = power_TDserial / (10^(EsN0 / 10));
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1i * randn(size(frame_TDserial)));
        frame_recieved = frame_TDserial + noise_msg;
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved, Nsym, Nframe* log2(M));
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        % QPSK解调
        frame_demod = QAMDemod(frame_recieved_FD, Nk, Nframe, M); %改为QAM解调
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
```

![在这里插入图片描述](https://img-blog.csdnimg.cn/direct/0acf84ed601f469695d81289cafb9f00.png#pic_center)


# 2 瑞利衰落信道之TDL信道

时延抽头延迟线(TDL，Tap Delay Line)信道是一种简化的多径模型，用于表示无线通道在特定延迟时间内的衰减特性。该模型通过一系列的时变复数增益，也称为"抽头"来表示信号经不同路径传播后的相位变化和衰减。在给定的TDL模型中，每个抽头都对应着一个特定的时延和相应的平均功率衰减。

信道的数学模型可以表示为：
$$
 y(t) = \sum_{i=1}^{N} h_i \cdot x(t - \tau_i) + n(t) 
$$


其中，
- \( y(t) \) 是接收信号，
- \( x(t) \) 是发送信号，
- \( h_i \) 是第 \( i \) 个抽头的复数增益，反映了信号在抽头 \( i \) 上的衰减和相位变化$h_i$ 取以 $ h_i = \sqrt{\text{PowerTDL}_i} \cdot \text{Rayleigh fading factor}$，并且 \( h \) 数组中的非零位置由TDL抽头的时延决定，则 $\text{Delay}[i]+1) = h_i $
- $ \tau_i$ 是第 i  个抽头的时间延迟，
- N是总抽头数，
- n(t) 是添加到通道输出的高斯白噪声。

在给定的信道模型中，PowerTDL_dB定义了各抽头的平均功率衰减，单位为分贝(dB)，Delay数组定义了各抽头相对于首个抽头的时延，单位为符号周期。这两个数组被用来计算PowerTDL，即各抽头的功率，并转换为线性单位。通过这些值，以及Rayleigh衰减模型，可以生成模拟多径效应的复数信道响应。

## 2.1 BPSK和QPSK



```python
%% 设置参数
clear;clc;
Nk = 128;           % 子载波个数
Nfft = 128;          % fft长度
Nframe = 6;         % 一帧中有几个OFDM符号
M = 1;              % 调制符号所含比特
SR = 250000;        % 符号速率
BR = SR .* M;       % 比特率
NGI = 32;           % 保护间隔长度
EbN0s = 0:2:20;      % 信噪比
Nsym = Nfft+NGI;    % 系统长度
bers = zeros(1,length(EbN0s));  % 误码率储存数组
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
    nloop = 10000;          % 发送多少帧
    n_biterror = 0;         % 错误的数据
    n_bitdata = 0;          % 一共发送了多少数据
    n_packeterror = 0;      % 有多少错误帧
    n_packetdata = 0;       % 发送了多少帧
    for ii = 1:nloop
%--------------------------发射端-------------------------------%
        % 生成一帧数据，串并转换，并QPSK，生成一帧
        
        frame_FDserial = randi([0 1], 1, Nk * Nframe * M);% 发送的是bit
        frame_FDparallel = reshape(frame_FDserial,Nk,Nframe*M);% 串并转换
        if M==1
            frame_mod = BPSKMod(frame_FDparallel, Nk, Nframe); % BPSK调制
        elseif M==2
            frame_mod = QPSKMod(frame_FDparallel, Nk, Nframe); % QPSK调制
            
        end
        % IFFT
        power_FT = sum(sum(abs(frame_mod).^2))/Nk/Nframe;  % 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_mod);         % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft);             % ifft
        % frame_ifft = ifft(frame_mod, Nfft);
        power_TD = sum(sum(abs(frame_ifft).^2))/Nk/Nframe; % 计算下IFFT前的能量，DT表示时域
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP");  % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI,1,Nsym*Nframe);
%--------------------------Channel-------------------------------%
        % 信号先经历衰落
        channel = Rayleigh_model(Nchannel, PowerTDL);
        h = zeros(1, Tau_maxTDL);
        h(Delay+1) = channel;
        frame_conv = conv(frame_TDserial, h);
        frame_fading = frame_conv(:,1:length(frame_TDserial));        % 看似是线性卷积，实际上由于CP变成了循环卷积
        % 添加高斯白噪声
        power_TDserial = sum(abs(frame_TDserial).^2)/Nk/Nframe;     
        EsN0 = EbN0 + 10*log10(M);                                  % 根据信噪比计算噪声能量，幅值，然后加在信号上
        N0 = power_TDserial .* 10.^(-EsN0/10);
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
        frame_recieved = frame_fading + noise_msg;
%--------------------------接收端-------------------------------%
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved,Nsym,Nframe);
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        % 信道均衡
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft));
        frame_equalization = frame_recieved_FD ./ repmat(H, 1, Nframe);
        
        if M==1
            % BPSK解调
            frame_demod = BPSKDemod(frame_equalization, Nk, Nframe);
        elseif M==2
            % QPSK解调
            frame_demod = QPSKDemod(frame_equalization, Nk, Nframe);
        end
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk*Nframe*M);
       
        % 计算error
        n_biterror_tmp = sum(abs(frame_output-frame_FDserial));
        n_bitdata_tmp = length(frame_FDserial);
        n_biterror = n_biterror + n_biterror_tmp;
        n_bitdata = n_bitdata + n_bitdata_tmp;
        if n_biterror_tmp ~= 0
            n_packeterror = n_packeterror + 1;
        end
        n_packetdata = n_packetdata + 1;
    end
    % 计算在当前信噪比下的误码率
    per = n_packeterror/n_packetdata;
    ber = n_biterror/n_bitdata;
    bers(kk)=ber;
    fprintf('%f\t%e\t%e\t%d\t\n',EbN0,ber,per,nloop);
end
%% 画图
semilogy(EbN0s,bers,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
% legend('理论曲线','实验曲线');
if M==1
    legend('BPSK实验曲线');
else
    legend('QPSK实验曲线');
end


function outs = QPSKMod(input_data,nk,nframe,m)
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*M)
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% m: 调制数
% 输出
% outs: (nk,nframe),输出a+bi
% out_coss:(nk,nframe)，输出实部
% out_sins:(nk,nframe)，输出虚部

if nargin < 4                   % 设置默认值
    m = 2;
end

outs = zeros(nk,nframe);    % 初始化
out_cos = zeros(nk,nframe);
out_sin = zeros(nk,nframe);

input_data_pn = 2*input_data - 1;     % 把0，1映射到-1，1
A = 1/sqrt(2);                        % 归一化幅值

out_cos((1:nk),(1:nframe)) = input_data_pn((1:nk), (1:2:2*nframe-1)) .* A;    % 每次使用两列
out_sin((1:nk),(1:nframe)) = input_data_pn((1:nk), ((2:2:2*nframe))) .* A;

outs = out_cos + out_sin * 1j;
end

function outputs = QPSKDemod(input_data,nk,nframe)
% 输入
% input_data: (Nk, Nframe), 一个频域的复数，会被拆开解调
% nk: 频域并联
% nframe: 一帧包含符号数
% 输出
% outputs：(Nk, 2*Nframe), 解调后，多出一倍，全是01
outputs = zeros(nk,2*nframe);
A = 1/sqrt(2); 
input_data = input_data ./ A;
outputs((1:nk),(1:2:2*nframe-1)) = real(input_data((1:nk),(1:nframe)))>=0;
outputs((1:nk),(2:2:2*nframe)) = imag(input_data((1:nk),(1:nframe)))>=0;

end

function outs = BPSKMod(input_data,nk,nframe,m)
% BPSK调制函数
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*M)
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% m: 调制数，默认值为1
% 输出
% outs: (nk,nframe),输出复数
A = 1/sqrt(2);    % 归一化幅值
outs = input_data * 2 - 1;  % 将二进制0/1映射为复数-1/1
outs = A * outs;  % 乘以归一化幅值
end

function outputs = BPSKDemod(input_data,nk,nframe)
% BPSK解调函数
% 输入
% input_data: (Nk, Nframe), 一个频域的复数
% nk: 频域并联
% nframe: 一帧包含符号数
% 输出
% outputs：(Nk, Nframe), 解调后的二进制数据
outputs = real(input_data) >= 0;  % 判断实部是否大于等于0，大于等于0表示1，小于0表示0
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
```
![在这里插入图片描述](https://img-blog.csdnimg.cn/direct/e18d46bcde0242849dd623fdcf122abe.png#pic_center)

## 2.2 QAM

```matlab
%% 设置参数
clear;clc;
Nk = 128; % 子载波个数
Nfft = 128; % fft长度
Nframe = 6; % 一帧中有几个OFDM符号
M = 4; % 调制符号所含比特（改为4：4QAM，16：16QAM,32：32:QAM）
SR = 250000; % 符号速率
BR = SR * log2(M); % 比特率（根据调制方式计算比特率）
NGI = 32; % 保护间隔长度
EbN0s = 0:1:12; % 信噪比
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
        power_TDserial = sum(abs(frame_TDserial).^2)/Nk/Nframe;     
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
```

![在这里插入图片描述](https://img-blog.csdnimg.cn/direct/729464c96bb24d0e811ff1004599b3be.png#pic_center)



