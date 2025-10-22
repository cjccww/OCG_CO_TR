function [y_out,nco_values]=symbolTimingPLL_Feedback_Fre( h,NFFT, x,  K1, K2, ...
     const, debug_s, debug_r)
% 频域方面  symbolPLL  算法

% 输入:
%   x       : 过采样接收信号, 复数矩阵 [L x nModes]
%   h       : 匹配滤波器冲激响应, 复数向量 [K x 1]
%   NFFT    : FFT长度
%   k1   : PI滤波器比例增益
%   k2   : PI滤波器积分增益
%
% 输出:
%   y_out      : 定时校正后的输出信号 [L x nModes]
%   nco_values : NCO相位值记录 [B x nModes]

%% 参数校验与默认值设置
if (nargin < 7)
    debug_s = 0; % 默认关闭静态调试绘图
end
if (nargin < 8)
    debug_r = 0; % 默认关闭实时调试示波器
end

if (debug_r)
    % 初始化图形
    figure;
    hFig = gcf;
    hAx = axes(hFig);
    hold on;
    grid on;
    box on;
    axis equal;
    % 设置坐标轴范围（类似原代码中的XLimits和YLimits）
    xlim(hAx, [-1.5 1.5]*max(real(const)));
    ylim(hAx, [-1.5 1.5]*max(imag(const)));
    xlabel('同相分量 (I)');
    ylabel('正交分量 (Q)');
    title('动态星座图显示');
    % 绘制星座点
    scatter(hAx, real(const), imag(const), 'x', 'filled', 'MarkerEdgeColor', 'r','LineWidth',2);
    % 初始化散点图对象用于动态更新
    hScatter = scatter(hAx, [], [], 'bo', 'filled');
    hScatter.XData=zeros(256,1);
    hScatter.YData=zeros(256,1);

end


if (debug_r)
  % 初始化图形
    figure;
    hFig = gcf;
    hFig.Name = 'Fractional Interval';  % 窗口标题
    hAx = axes(hFig);
    hold on;
    grid on;  % 显示网格
    box on;
    xlabel('Sample Index');
    ylabel('Amplitude');
    title('Fractional Interval');
    % 用于存储最近的波形数据（对应BufferLength）
    bufferLength = 1e2;
    waveformBuffer = zeros(1, bufferLength);
    currentIndex = 0;
    % 设置坐标轴范围
    xlim(hAx, [1 bufferLength]);  % 时间跨度1e4个点
    ylim(hAx, [-1 1]);   % Y轴范围[-1,1]

    % 初始化波形对象 - 使用初始点避免空对象
    hLine = plot(hAx, 1, 0, 'b');  % 先绘制第一个点
    
end

%% 数据预处理；参数设定

% 调制阶数提取，输入信号格式化，输入信号选择
M = numel(const); % 根据星座图获取调制阶数

% 将输入信号转为列向量（确保后续处理一致性）
if (size(x, 1) == 1)
    x = x(:);
end
L = size(x, 1); % 信号长度
K = length(h);  % 滤波器长度
delay = floor((K - 1) / 2); % 滤波器群延迟

if NFFT < K
    error('NFFT deve ser maior ou igual ao comprimento do filtro');
end

% --- 预计算 ---
overlap = K - 1; % 重叠保留法的重叠点数
f1 = 0:1/NFFT:0.5-1/NFFT; % 生成频率向量 [0, 1/NFFT, ..., 0.5*(NFFT-1)/NFFT]
f2=-0.5:1/NFFT:0-1/NFFT;  % [-0.5*(NFFT-1)/NFFT,1/NFFT, ...,0]
f=[f1,f2];

% 计算总块数 B
B = ceil((L + overlap) / (NFFT - K + 1)); % NFFT - overlap = NFFT - K + 1

% ------------------ 初始化缓冲区 ---------------
% 匹配滤波器填充
h=[h;zeros(NFFT-K, 1)];
% N: 信号块缓冲区
% 我们只处理NFFT长度的块
N = zeros(NFFT+delay, 1); % 当前信号块

% 计算输出信号长度
output_length = B * (length(N) - overlap);
y_out = zeros(output_length, 1); % 预分配输出
nco_values = zeros(length(x), 1); % 存储NCO输出相位
%%  TimingRecovery feedback
% 在前面填充 overlap 个零，在后面填充 NFFT 个零
signalPad = [zeros(overlap, 1); x; zeros(NFFT, 1)];

% --- 初始化DPLL状态 ---
out_nco = 0;    % NCO当前输出相位
intPart = 0;    % PI滤波器积分部分

% --- 块处理循环 ---
for indB = 1:B

    % --- 1. 提取信号块 ---
    step_start = (indB - 1) * (NFFT - overlap) + 1; % MATLAB索引从1开始
    end_step = step_start + NFFT - 1;

    signal = signalPad(step_start:end_step);
    % --- 2. 频域处理 ---
    outBlockFFT = fft(signal); % FFT
    H=fft(h);% 信道的频域函数

    % --- 3. 频域相位旋转 (插值) ---
    % 使用 exp(-1j*2*pi*f*out_nco), f 是频率向量
    % out_nco 是归一化的相位 (单位: 符号周期)
    phase_rotation = exp(-2j*pi*f'*out_nco); % f' 是列向量
    Y = outBlockFFT.* H .* phase_rotation;

    % --- 4. Godard TED ---
    method = 'Modity'; % 'Modity' 'Conventional'
    ted(indB) = godardTED(Y, NFFT,method);

    % --- 5. PI 环路滤波器 ---
    intPart = intPart + K1 * ted(indB); % 积分部分更新
    propPart = K2 * ted(indB);          % 比例部分
    loopFilterOut(indB) = propPart + intPart;
    
    % --- 6. NCO 更新 ---
    out_nco = out_nco - loopFilterOut(indB); % 更新相位

    % --- 7. NCO 相位归一化 (模1) ---
    % (处理NCO相位溢出/欠溢，管理采样索引)
    while out_nco > 1
        % 如果NCO相位超前超过1个采样间隔
        out_nco = out_nco - 1;
    end
    while out_nco < -1
        % # 如果NCO相位滞后超过1个采样间隔
        out_nco = out_nco + 1;
    end

    % --- 8. IFFT 回到时域 ---
    outBlockIFFT = ifft(Y);

    % --- 9. 保存有效输出 (重叠保留法) ---
    % 丢弃前面 overlap 个点
    valid_start = overlap + 1;
    valid_end = NFFT;
    valid_data = outBlockIFFT(valid_start:valid_end);

    % 计算在y_out中的写入位置
    y_out(step_start:step_start+(NFFT-overlap)-1)=valid_data;
    % --- 10. 记录NCO值 ---
    nco_values(step_start:step_start+NFFT-1) = out_nco;


    if (debug_r)
         %  更新星座图（固定长度数据点）
            hScatter.XData = real(y_out(step_start:step_start+(NFFT-overlap)-1));
            hScatter.YData = imag(y_out(step_start:step_start+(NFFT-overlap)-1));


            % 更新实时μ值曲线
            currentIndex = currentIndex + 1;
            if currentIndex > bufferLength
                % 当缓冲区满时，移除最旧的数据，添加新数据
                waveformBuffer = [waveformBuffer(2:end), nco_values(step_start)];
                xlim(hAx, [currentIndex-bufferLength+1 currentIndex]);
            else
                % 缓冲区未满时，直接添加新数据
                waveformBuffer(currentIndex) = nco_values(step_start);
            end
            if currentIndex < bufferLength
                % 确定当前显示的范围（最近1e4个点）
                displayStart = max(1, currentIndex - bufferLength + 1);
                displayData = waveformBuffer(displayStart:currentIndex);
                displayX = displayStart:currentIndex;
            else
                displayData = waveformBuffer;
                displayX=currentIndex -bufferLength + 1:currentIndex;
            end

            % 使用下标赋值方式更新数据
            hLine.XData(1:length(displayX)) = displayX;
            hLine.YData(1:length(displayX)) = displayData;

            % 刷新图形
            drawnow limitrate;  % 限制刷新速率，提高性能
            pause(0.2)
    end
end % for indB
% --- 11. 最终输出对齐 (补偿滤波器延迟) ---
% 由于滤波器h引入了delay个采样点的延迟，需要截取
% 处理流程中，重叠保留法本身也引入了延迟。
% 通常，最终输出需要截取从 'delay' 开始的 L 个点。
if size(y_out, 1) >= delay + L
    y_out = y_out(delay+1:delay+L);
else
    % 如果输出不够长，可能需要填充零或报错
    error('输出信号长度不足，无法对齐！');
end

%% 参数可视化

if (debug_s)
    figure
    plot(ted)
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')

    figure
    plot(loopFilterOut)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(nco_values, '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
end
