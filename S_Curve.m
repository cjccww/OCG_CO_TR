% CO system Train
clc;close all;clear;
addpath('Phase_Sync\')
addpath('Dsp\')
addpath('Fncs\')

%% 系统初始化
% 信号生成
coherentGeneration;
% 上采样率
SPS = Tx.TxPHY.sps;
% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);

% 发射信号 双偏振
[signalDualPol,qamDualPol]  =  Tx.dataOutput();

% DSP参考信号
label=qamDualPol;

% 根升余弦发射滤波器（Tx）
htx = rcosdesign(Tx.Nr.psfRollOff, RcDelay, SPS);
% 匹配滤波器（Rx，时域反转共轭）
mf = conj(fliplr(htx));
% 组合升余弦脉冲（Tx+Rx级联响应）
r_p = conv(htx, mf);
% 数据上采样
symbolsUp = upsample(qamDualPol, SPS);
% 匹配滤波器输出（含成形滤波）
for indMode = 1:Tx.TxPHY.Nmodes
    mfOutput1(:,indMode) = conv(r_p, symbolsUp(:,indMode));
end

% 数据性能起始位置
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
%% Tx and Rx Lo
% phase noise TX
lw   =  100e3;
RIN  =  0;
Plo_dBm  =  10;
phase_Noise  =  Tx.phaseNoiseMode(signalDualPol,lw);
% Rx LO
sigLO   = Tx.laserMode(signalDualPol,lw,RIN,Plo_dBm);
% 添加频偏
FO      = 0;                 % frequency offset
sigLO_Offset = Tx.addFrequencyOffset(sigLO,FO); % add frequency offset

%% 器件频响建立

% obj.Implementation.responType='Bessel';
% f3dB=20e9;
% order=5;
% verbose=0;
% [filt,H]=Rx.createFrequencyResponse(freq,order,f3dB,verbose);
% out1 = firFilter(filt.h, txSig); % 第一种滤波器工作方式,时域
% out2=filter_work(txSig,H); % 第二种滤波器工作方式，频域

%% 双偏发射机创建
Amp=0.5; % EA放大
for indMode = 1:Tx.TxPHY.Nmodes

    % X,Y偏振信号传输
    if indMode==1
        sigTxCh = iqm(phase_Noise, Amp*signalDualPol(:,indMode).',paramIQ);
    else
        sigTxCh = DP_iqm(phase_Noise, Amp*signalDualPol(:,indMode).',paramIQ);
    end
    % 每个通道的输出模式
    fprintf(' mode %d\t power: %.2f dBm\n', indMode, 10 * log10(((10 .^ (Pout_dBm(indMode) / 10)) * 1e-3/Tx.TxPHY.Nmodes) / 1e-3));
    % 设置功率
    sigTxCh=Tx.setSignalPower(sigTxCh,channelPowerType,Pout_dBm(inMode));
    power=signalpower(sigTxCh);
    fprintf(' mode %d optical signal power: %.2f dBm\n',indMode, 10 * log10(power / 1e-3));
    % 装载信号
    sigDualPolTx(:,indMode)=sigTxCh;

end

%% 添加信道延迟
% timeOffset = 25;       % 信道延迟（采样点数）
% % 延迟模块（模拟信道延迟）
% DELAY = dsp.Delay(timeOffset);
% % 信道延迟
% delaySig = step(DELAY, sigDualPolTx);    % 添加固定延迟;直接添加零

% % 延迟
% delay=[38*1e-12,35*1e-12];
% for indMode=1:Tx.TxPHY.Nmodes
%     delaySig(:,indMode)=Rx.addSkew(sigDualPolTx(:,indMode),delay(indMode));
% end

%% 线性信道
CD = 0;
paramFiber = struct();
paramFiber.L = 100  ;
paramFiber.alpha = 0 ;
paramFiber.Fs = fs  ;
paramFiber.Fc = Fc  ;
paramFiber.D = CD;
for indMode = 1:Tx.TxPHY.Nmodes
    sigRxo(:,indMode)= linearChannel(sigDualPolTx(:,indMode), paramFiber);
end

%% PMD 及 RSOP 
% method 1 
[f, t_pmd]=freq_time_set(size(sigRxo,1),fs);
w = 2*pi * f.';
% PDL损耗
meanPDLdB=1;             %引入PDL相关损耗
gama=10^(-meanPDLdB/20);    %注意在这里加入了DGD损耗
% DGD设置—— 以一个符号周期为标尺尺度
DGD_symbl = 1/fb;
meanDGD = DGD_symbl*0.4; % 平均的DGD设置为0.4个符号周期
PP1= exp(1j*w*meanDGD/2);  % 快轴相移，对角相，引入DGD
PP2 = 0;
% RSOP旋转
rsop_theta = pi/8;
% RSOP 旋转
mat_rot = [cos(rsop_theta),-sin(rsop_theta);...
    sin(rsop_theta),cos(rsop_theta)];

% 将X和Y偏振信号转换到频域
sig_fft_x = fft(sigRxo(:,1));
sig_fft_y = fft(sigRxo(:,2));
% 将零频率分量移到频谱中心
sig_fft_x = fftshift(sig_fft_x);
sig_fft_y = fftshift(sig_fft_y);

% PMD添加
uux = gama.*PP1.* sig_fft_x + PP2.* sig_fft_y;
uuy = -conj(PP2).* sig_fft_x + conj(PP1).*sig_fft_y;
% 转换为时域
sig_x=ifft(fftshift(uux));            %%这里对应下文的re_ex
sig_y=ifft(fftshift(uuy));

% RSOP 应用
sig_x_out = mat_rot(1, 1) .* sig_fft_x + mat_rot(1, 2).* sig_fft_y;
sig_y_out = mat_rot(2, 1).* sig_fft_x + mat_rot(2, 2).* sig_fft_y;
% 返回处理后的信号
sigRxo_out(:,1) = sig_x_out;
sigRxo_out(:,2) = sig_y_out;

%% 使用先变化正交基，添加PMD，再将正交基转回
% method 2  
% 只使用一段
pmd_trunk_num = 1 ;
% 旋转角度设置
psp_manual = 1; % 1代表随机数 0代表分布
if psp_manual
    rsop_theta = pi/4;
    psp_theta = rsop_theta* ones(pmd_trunk_num, 1);
    psp_phi = rsop_theta * ones(pmd_trunk_num, 1);
else
    % 均匀分布在[-π, π)区间
    psp_theta = rand(pmd_trunk_num, 1) * 2 * pi - pi;
    % 在Poincare球面上的均匀分布
    psp_phi = 0.5 * asin(rand(pmd_trunk_num, 1) * 2 - 1);
end

% DGD设置—— 以一个符号周期为标尺尺度
DGD_symbl = 1/fb;
meanDGD = DGD_symbl*0.4; % 平均的DGD设置为0.4个符号周期
sigma =  0.001; % 方差
DGD_rand=sqrt(3*pi/(8*pmd_trunk_num))*(1+sigma*rand(pmd_trunk_num,1)).*meanDGD;
pmd_beta = w .* DGD_rand/2;

% 定义Pauli矩阵
sig0 = eye(2);  % 单位矩阵
sig2 = [0, 1; 1, 0];  % Pauli σ2矩阵
sig3i = [0, 1; -1, 0];  % Pauli σ3矩阵的虚数倍

% 构建旋转矩阵的一部分（θ旋转）
mat_theta = cos(psp_theta) * sig0 ...
    - sin(psp_theta) * sig3i ;  % 正交矩阵
% 构建旋转矩阵的另一部分（φ旋转）
mat_epsilon = cos(psp_phi) * sig0 ...
    + 1j * sin(psp_phi) * sig2;  % 正交矩阵

% 组合成完整的旋转矩阵：用于将信号变换到主态(PSP)基
mat_rot = mat_theta * mat_epsilon;  % 得到RSOP矩阵
mat_rot_conj = conj(mat_rot);  % 旋转矩阵的共轭

% 注意：设A=[Ax;Ay]为电场，则matR*D*matR'*A是线性PMD步骤，其中D是DGD矩阵

% 1> 通过乘以旋转矩阵R的共轭，将信号变换到PSP基上
uux = mat_rot_conj(1, 1) * sig_fft_x + mat_rot_conj(2, 1) * sig_fft_y;
uuy = mat_rot_conj(1, 2) * sig_fft_x + mat_rot_conj(2, 2) * sig_fft_y;

% 对两个主态分量应用不同的相位延迟
uux = exp(-1j * (pmd_beta)) .* uux;
uuy = exp(-1j * (-pmd_beta)) .* uuy;
% 3> 通过乘以旋转矩阵R，将信号变换回原始基
sig_fft_x_method2 = mat_rot(1, 1) * uux + mat_rot(1, 2) * uuy;
sig_fft_y_method2 = mat_rot(2, 1) * uux + mat_rot(2, 2) * uuy;

% 将频谱移回标准格式 ,并转换为时域
sig_x_method2 =ifft(ifftshift(sig_fft_x_method2)) ;
sig_y_method2 = ifft(ifftshift(sig_fft_y_method2));

% 返回处理后的信号
sigRxo_out_method2(:,1) = sig_x_method2;
sigRxo_out_method2(:,2) = sig_y_method2;


%% 信号传输
% sigRxo=Rx.signalTran(sigDualPolTx,param);
% if Tx.TxPHY.Nmodes==2
%     fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
%     fprintf('total ssfm signal power model 2: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,2))/ 1e-3));
% else
%     fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
% end
% % 频谱参考
% mon_ESA(sigRxo,fs);

%% 相干接收
theta=0;
sig_RxE=Rx.coherentReceive(sigRxo,sigLO_Offset,theta,paramPD);

%% 应用集成的时钟偏移函数 及 timing jitter
% ppm=-20;
% jitter_rms=0; %1e-10
% txResamp1=Rx.addSamplingClockOffset(txSig,ppm,jitter_rms);

%% 匹配滤波
matchOut1=Rx.matchFiltering(sig_RxE);
%% 加噪
Eb_N0_dB = 10;
for i  = 1:50
    % 选择一个偏振态
    matchOut=Rx.addNoiseEbN0(matchOut1(:,1).',Eb_N0_dB);
    matchOut=matchOut.';
    %% S曲线 延迟设置
    tau = 0;  % 假设真实定时偏移为0（基准点）
    tauEstVec = -(SPS/2):(SPS/2);      % 接收机估计的偏移范围：-L/2 ~ L/2 采样点
    tauErrVec = tau - tauEstVec;   % 实际误差 = 真实偏移 - 估计偏移

    % 偏移归一化
    normTauE = tauErrVec/ SPS;
    % 理想符号采样点索引
    idealStrobeIdx =(RcDelay * SPS)+ (1:SPS:nSymbols*SPS);
    idealStrobeIdx_2sps = (RcDelay * SPS)+ (1:SPS/2:nSymbols*SPS);
    idealStrobeIdx_4sps =  (RcDelay * SPS)+(1:SPS/4:nSymbols*SPS);
    % 创建参数
    average = zeros(length(tauEstVec),1);

    for ind=1:length(tauEstVec)
        tauEst = tauEstVec(ind);      % 接收机假设的定时偏移（采样点）
        tauErr = tau + tauEst;      % 实际定时误差

        % 实际采样点索引（引入估计偏移）
        strobeIdx_1sps = idealStrobeIdx + tauEst;
        strobeIdx_2sps=idealStrobeIdx_2sps+ tauEst;
        strobeIdx_4sps=idealStrobeIdx_4sps+ tauEst;

        %%
        TED = 'lee';
        switch (TED)

            case 'gardner'
                zcStrobeIdx = strobeIdx_1sps(2:end) - SPS/2;   % 符号间过零点位置
                prevStrobeIdx = strobeIdx_1sps(1:end-1);     % 前符号位置
                currStrobeIdx = strobeIdx_1sps(2:end);       % 当前符号位置
                % 过零点采样 × (前符号 - 当前符号)
                err = real( conj(matchOut(zcStrobeIdx)) .* ...
                    (matchOut(prevStrobeIdx) - matchOut(currStrobeIdx))); % 也是符合公式
            case 'mm' % Mueller-Muller TED
                prevStrobeIdx = strobeIdx_1sps(1:end-1);     % 前符号位置
                currStrobeIdx = strobeIdx_1sps(2:end);       % 当前符号位置
                % 前符号 × 当前采样 - 当前符号 × 前采样
                err = qamDualPol(1:end-1,1) .* matchOut(currStrobeIdx) - ...
                    qamDualPol(2:end,1)   .* matchOut(prevStrobeIdx);

            case 'gardnerNyquist'

                zcStrobeIdx = strobeIdx_1sps(2:end) - SPS/2;   % 符号间过零点位置
                prevStrobeIdx = strobeIdx_1sps(1:end-1);     % 前符号位置
                currStrobeIdx = strobeIdx_1sps(2:end);       % 当前符号位置
                % 功率型 Gardner-Nyquist
                err = abs(matchOut(zcStrobeIdx)).^2 .* (abs(matchOut(currStrobeIdx)).^2 - abs(matchOut(prevStrobeIdx)).^2);

            case 'godard'
                NFFT = 2048;
                sigRxOffset = matchOut(strobeIdx_2sps);
                % 零填充到 FFT 长度的整数倍
                lenPad = NFFT * ceil(numel(sigRxOffset)/NFFT);
                lpad = [sigRxOffset(:); zeros(lenPad - numel(sigRxOffset), 1)];

                numBlocks = lenPad / NFFT;
                err = zeros(1, numBlocks);

                for blk = 1:numBlocks
                    startIdx = (blk-1)*NFFT + 1;
                    endIdx   = blk*NFFT;
                    blkSig = lpad(startIdx:endIdx);

                    blkFFT = fft(blkSig);
                    err(blk) = godardTED(blkFFT,NFFT,'Modity');
                end
            case 'sln'
                % 块的长度
                sigRxOffset=matchOut(strobeIdx_4sps);
                NFFT = length(sigRxOffset);
                % --- 步骤2: 计算块的数量 ---
                numBlocks = floor(length(sigRxOffset) / NFFT); % 使用floor确保是整数
                % --- 步骤3: 初始化误差存储 ---
                err = zeros(numBlocks, 1); % 预分配列向量，提高效率
                for blockIdx=1:numBlocks
                    px=sigRxOffset((blockIdx-1)*NFFT+1:blockIdx*NFFT);
                    N = length(px);
                    k = 1:N;
                    ex = exp(-1j.*(k-1).*pi./2);
                    s = sum( abs(px).^2 .* ex.' );
                    temp = 1/2/pi*angle(s);
                    err(blockIdx)=temp;
                end

            case 'lee'
                sigRxOffset = matchOut(strobeIdx_2sps);
                NFFT = length(sigRxOffset);
                g_bias=1; % 偏执项
                % --- 步骤2: 计算块的数量 ---
                numBlocks = floor(length(sigRxOffset) / NFFT); % 使用floor确保是整数
                % --- 步骤3: 初始化误差存储 ---
                err = zeros(numBlocks, 1); % 预分配列向量，提高效率

                for blockIdx=1:numBlocks
                    px=sigRxOffset((blockIdx-1)*NFFT+1:blockIdx*NFFT);
                    L1 = length(px);
                    n = 1:L1;
                    % cosine part
                    ex1 = (-1).^(n-1);          % ex1 = exp(-1j.*(ii-1).*pi);
                    sum_1 = sum( abs(px).^2 .* ex1.' );

                    % sine part
                    ex2 = 1j * (-1).^(n-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
                    xh = px(2:end);
                    xx = px(1:end-1);
                    ex2 = ex2(1:end-1);
                    sum_2 = sum( real(conj(xx).*xh) .* ex2 .');

                    ex3 = exp(-1j.*(n-2).*pi);
                    ex3 = ex3(1:end-1);
                    sum_3 = ( sum ( (xx+1j*xh).* (conj(xx)+1j*conj(xh)) .* ex3.') );

                    % with biasing
                    s = g_bias*sum_1 + sum_2;
                    %s = conj(s);

                    temp = 1/2/pi * angle(s);
                    % 近似使用
                    temp_1 =1/2/pi * angle(sum_3);
                    err(blockIdx)=temp;
                end


        end
        average(ind) = mean(err);
    end

    Total_average(:,i)=average;

end

%% 结果可视化
figure;
plot(normTauE,Total_average)
% 搜寻所需点
req_Tao = getOSNRPenalty(normTauE(end-5:-1:5), Total_average(5:end-5,:), 0);
% histogram(req_Tao)

jitter_dB=20*log10(sqrt(var(req_Tao)));
fprintf('Timing Jitter: %.2f dBm\n',jitter_dB);


% 绘制统计直方图
scurve_v = mean(Total_average, 2);
diff_matrix = Total_average - scurve_v;
figure;
histogram(real(diff_matrix), 10, 'Normalization', 'count');
title('频数统计');
figure;
histogram(real(diff_matrix), 10, 'Normalization', 'pdf');
title('概率密度统计');
%%
function errors=godardTED(blk,NFFT,method)
% 频域TED实现
switch (method)

    case 'Modity'
        % --- 步骤4.3: 计算Godard TED误差 ---
        error_val = sum(imag( blk(1:floor(NFFT/2)) .* conj( blk(floor(NFFT/2)+1:end) ) ) );
    case 'Conventional'
        % 两种误差计算方法
        error_val= (1/(2*pi)) * angle( sum( blk(1:NFFT/2) .* conj( blk(NFFT/2+1:end) ) ) );
end
% 存储误差
errors = error_val/NFFT;

end % function