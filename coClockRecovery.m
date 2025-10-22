% CO system Train
clc;close all;clear;

%% 信号生成
coherentGenerationClock;
% 添加路径
addpath('Plot\')
addpath('Dsp\')
addpath('Phase_Sync\')
addpath('Fncs\')


%% 调试配置
debug_tl_static  = 0; % 静态调试标志（1=显示最终星座图）
debug_tl_runtime = 0; % 运行时调试标志（1=显示同步过程示波器）
debug_tl_runtime_fast = 1; % 快速运行调试状态
%% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% 发射信号 
[txSig,qamSig]  =  Tx.dataOutput();

% 参考信号
label=qamSig;
% 数据性能起始位置
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
%% 定时恢复环路参数
Bn_Ts    = 0.01;       % 环路带宽×符号周期（归一化噪声带宽）
eta      = 1;          % 环路阻尼系数（控制收敛速度）

% 信号处理参数
timeOffset = 0;       % 信道延迟（采样点数）
fsOffsetPpm =100;       % 采样时钟频偏（ppm单位）
EsN0     = 20;         % 信噪比（Es/N0 dB）
Ex       = 1;          % 符号平均能量
TED      = 'GTED';    % 定时误差检测器类型（关键算法选择） 'ELTED', 'ZCTED', 'GTED', 'MMTED','MLTED'% 五种算法
intpl    = 3;          % 插值方法：0=多相，1=线性，2=二次，3=三次
forceZc  = 0;          % 强制零交符号模式（调试自噪声）
% 延迟模块（模拟信道延迟）
DELAY = dsp.Delay(timeOffset);

% 参考信号

M=Tx.TxPHY.M;
const = qammod(0:M-1,M );          % QAM星座图
Ksym = modnorm(const, 'avpow', Ex);   % 能量归一化因子
const = Ksym * const;                 % 调整后的参考星座
% 调制误差率(MER)测量模块
mer = comm.MER;
mer.ReferenceSignalSource = 'Estimated from reference constellation';
mer.ReferenceConstellation = const;   % 设置参考星座

%% 定时恢复环路常数计算
% 计算定时误差检测器增益
Kp = calcTedKp(TED, Tx.Nr.psfRollOff);         % 自定义函数计算TED增益
% Kp = calcTedKp('Lee',Tx.Nr.psfRollOff,'simulated');
% 调整增益（考虑符号能量）
K  = 1;           % 假设信道增益为1（实际系统需AGC）
Kp = K * Ex * Kp; % 调整后的TED增益

% 定时环路控制参数
K0 = -1;          % 计数器增益（类似DDS的灵敏度）
[ K1, K2 ] = piLoopConstants(Kp, K0, eta, Bn_Ts, Tx.TxPHY.sps); % PI控制器参数

% 打印环路参数
fprintf("环路常数:\n");
fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);
%% ClockRecover paramer
% 参数初始化
clockRecovery = DspSyncDecoding( ...
    fs,...         % 接收信号的采样率
    fb, ...        % 接收信号的波特率
    16,...         % 接收信号的格式
    fs/fb, ...     % 上采样率
    [],...         % 时钟信号的上采样率
    skip, ...      % 误码计算起始位置
    label,...      % 参考信号
    "MLTED");             
clockRecovery.Implementation.eta       = eta;              % 环路阻尼因子（稳定性控制）
clockRecovery.Implementation.Bn_Ts     = Bn_Ts;           % 环路带宽 × 符号周期（控制同步速度）
clockRecovery.Implementation.Ex =1;
clockRecovery.Implementation.intpl=2;          % 插值方法：0=多相，1=线性，2=二次，3=三次
%% Train
% 模拟采样时钟偏移（通过重采样）【可以进行优化，使用插值函数】
fsRatio = 1 + (fsOffsetPpm * 1e-6); % 计算频率比例（Rx/Tx）
tol = 1e-9;                          % 重采样容差
[P, Q] = rat(fsRatio, tol);          % 转换为有理分数
txResamp = resample(txSig, P, Q);    % 执行重采样


% % 信号延迟周期采样点
% delay=ta/2;
% % 集成的延迟函数
% % delaySig=Rx.addSkew(txResamp1,delay);
% txSig=Rx.complex_signal_delay(txSig,delay);

% 应用集成的时钟偏移函数
ppm=-5;
jitter_rms=0; %1e-10
txResamp1=Rx.addSamplingClockOffset(txSig,ppm,jitter_rms);

% % 信号延迟周期采样点
% delay=ta/2;
% % 集成的延迟函数
% % delaySig=Rx.addSkew(txResamp1,delay);
% delaySig=Rx.complex_signal_delay(txResamp1,delay);

delaySig=txResamp1;
% delaySig = step(DELAY, txResamp1);    % 添加固定延迟;直接添加零

% 添加高斯白噪声（AWGN）
txSigPower = 1 / sqrt(Tx.TxPHY.sps);            % 计算信号功率
rxSeq = awgn(delaySig, EsN0, txSigPower); % 添加带限噪声

% match 
mfOut=Rx.matchFiltering(rxSeq);

%% 定时恢复处理
% 无定时恢复的直接抽取
rxNoSync = downsample(mfOut, Tx.TxPHY.sps);     % 简单抽取（存在时偏）

% 理想定时恢复（已知延迟）
rxPerfectSync = downsample(mfOut, Tx.TxPHY.sps, timeOffset); % 完美同步

%定时恢复算法
rollOff=Tx.Nr.psfRollOff;
rcDelay=Tx.Nr.psfLength ;
sps=Tx.TxPHY.sps;
type = 'feedback';
% type = 'feedford';
% TED='SLNTED';  % 'SLNTED''LeeTED'
% [ rxSync1 ] = symbolTimingSync(TED, intpl, sps, rxSeq, mfOut, K1, K2, ...
%     const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime);

% 适用于SLN参数
% K2 = 1e-3;
% K1 = -4e-3; % 4e5
[ rxSync2 ]=symbolTimingPLL(TED, intpl, sps, rxSeq, mfOut, K1, K2, ...
    const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime,debug_tl_runtime_fast,type);
y2=downsample(rxSync2,sps);


% 频域算法
h=ones(1,1);
NFFT=1024;

[y_out,nco_values]=symbolTimingPLL_Fre( h,NFFT, mfOut,  K1, K2, ...
     const, 0, 0);
y1=downsample(y_out,sps);

%% 性能评估


% 计算并显示MER结果
fprintf("\n测量MER结果:\n")
fprintf("无定时恢复: %.2f dB\n", mer(rxNoSync(skip:end)))
fprintf("理想定时恢复: %.2f dB\n", mer(rxPerfectSync(skip:end)))
fprintf("频域时钟恢复算法 %s算法: %.2f dB\n", TED, mer(y1(skip:end)))
fprintf("时钟恢复算法 %s算法: %.2f dB\n", TED, mer(y2(skip:end)))

%% 调试绘图（需要手动开启debug标志）

if (debug_tl_static)
     scatterplot(rxNoSync(skip:end));      % 未同步星座
    title('无定时恢复');

     scatterplot(rxPerfectSync(skip:end)); % 理想同步
    title('理想定时恢复');

     scatterplot(rxSync2(skip:end));       % 
    title(sprintf('时钟恢复算法 %s 算法', TED));
end