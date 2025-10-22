function Filtered = OptiFilter(Signal, SymbolRate, Co_band, Nt, NSymbol, Ftype, Bb)
% OptiFilter: 光/电信号频域滤波器（高斯、超高斯、6阶Butterworth、5阶Bessel）
%
% 输入:
%   Signal      : 输入信号向量（列向量，复包络或实信号）
%   SymbolRate  : 符号速率 [Hz]
%   Co_band     : 相对带宽系数，-3dB带宽 = SymbolRate * Co_band
%   Nt          : 每符号采样点数（过采样率）
%   NSymbol     : 符号总数
%   Ftype       : 滤波器类型（字符串）
%                 'gaussian'        : 高斯滤波器
%                 'supergaussian4'  : 4阶超高斯
%                 'butter6'         : 6阶Butterworth
%                 'bessel5'         : 5阶Bessel
%   Bb          : Bessel滤波器归一化系数（仅对 'bessel5' 有效）
%
% 输出:
%   Filtered    : 滤波后的时域信号，长度与Signal一致
%
% 注意:
%   - 本函数使用FFT/IFFT在频域完成滤波，适用于大长度信号
%   - 所有滤波器响应均为解析式，无需额外设计
%   - 支持复信号输入，适用于相干光通信系统
%
% 示例:
%   Filtered = OptiFilter(Signal, 32e9, 0.75, 32, 1024, 'gaussian', 0);
% 参数示例
% Bb = 0.28;
%% --------------------- 参数检查 ---------------------
if nargin < 6
    error('输入参数不足，请检查输入。');
end

if ~ischar(Ftype)
    error('Ftype必须为字符串。');
end

Ftype = lower(Ftype);  % 统一小写，避免大小写错误

%% --------------------- 滤波器系数定义 ---------------------
% 6阶Butterworth系数（已预计算）
bu1 = 3.86370330515627315;
bu2 = 7.4641016151377546;
bu3 = 9.1416201726856413;
bu4 = bu2;
bu5 = bu1;

% 5阶Bessel系数（分子分母）
be0 = 945;
be1 = 945;
be2 = 420;
be3 = 105;
be4 = 15;

%% --------------------- 频域参数设置 ---------------------
N = Nt * NSymbol;               % 总采样点数
T = 1 / SymbolRate;             % 符号周期
Ts = T / Nt;                    % 采样周期
Fs = 1 / Ts;                    % 采样率
DeltaF = Fs / N;                % 频率分辨率

bw = SymbolRate * Co_band;      % -3dB带宽（Hz）
f = SymbolRate / NSymbol * [-N/2 : N/2 - 1]';  % 频率轴（列向量）
x = f / bw;                     % 归一化频率（相对于带宽）

%% --------------------- 初始化输出 ---------------------
Filtered = zeros(N, 1);         % 预分配输出向量

%% --------------------- 构造频域响应H(f) ---------------------
switch Ftype
    case 'gaussian'
        % 高斯滤波器：exp(-0.5 * ln(2) * x^2)
        Hf = exp(-0.5 * log(2) * x.^2);

    case 'supergaussian4'
        % 4阶超高斯滤波器：exp(-0.5 * ln(2) * x^(2*4))
        Hf = exp(-0.5 * log(2) * x.^(2*4));

    case 'butter6'
        % 6阶Butterworth滤波器（复响应）
        x2 = x.^2;
        x3 = x2 .* x;
        x4 = x3 .* x;
        x5 = x4 .* x;
        x6 = x5 .* x;
        Hf = 1 ./ (1 - bu2*x2 + bu4*x4 - x6 + 1i*(bu1*x - bu3*x3 + bu5*x5));

    case 'bessel5'
        % 5阶Bessel滤波器（复响应）
        om  = 2 * pi * x * Bb;
        om2 = om.^2;
        om3 = om2 .* om;
        om4 = om3 .* om;
        om5 = om4 .* om;

        pre = be0 - be2*om2 + be4*om4;  % 实部
        pim = be1*om - be3*om3 + om5;   % 虚部
        Hf = be0 ./ (pre + 1i*pim);

    otherwise
        error('不支持的滤波器类型: %s', Ftype);
end

%% --------------------- 频域滤波 ---------------------
Hf = fftshift(Hf);                          % 将0频移至中心，匹配FFT格式
Freq_Signal = fft(Signal);                  % FFT变换
Filtered_freq = Freq_Signal .* Hf;          % 频域相乘
Filtered = ifft(Filtered_freq);             % IFFT回时域

end