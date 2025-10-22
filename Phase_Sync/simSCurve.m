function [ normTauE, g ] = simSCurve(TED, rollOff, rcDelay, nSymbols)
%% 功能：通过仿真生成定时误差检测器（TED）的S曲线
% 输入：
%   TED      - 定时误差检测器类型（MLTED/ELTED/ZCTED/GTED/MMTED）
%   rollOff  - 升余弦滚降因子
%   rcDelay  - 升余弦滤波器时延（符号数，默认10）
%   nSymbols - 仿真符号数量（默认1e4）
% 输出：
%   normTauE - 归一化定时误差向量（范围[-0.5, 0.5]，单位符号周期）
%   g        - S曲线值向量（TED输出随定时误差的变化）

%% 参数初始化
if nargin < 3
    rcDelay = 10;    % 默认升余弦滤波器时延为10个符号
end
if nargin < 4
    nSymbols = 1e4;  % 默认仿真1万符号
end

%% 滤波器设计
L = 1e3;  % 过采样因子（仅用于高精度仿真，与实际通信系统无关）

% 根升余弦发射滤波器（Tx）
htx = rcosdesign(rollOff, rcDelay, L);

% 匹配滤波器（Rx，时域反转共轭）
mf = conj(fliplr(htx));

% 组合升余弦脉冲（Tx+Rx级联响应）
r_p = conv(htx, mf);

% 导数匹配滤波器（用于MLTED）
dmf = derivativeMf(mf, L);

%% 信号生成与处理
% 生成2-PAM随机符号（数据辅助模式）
data = randi(2, 1, nSymbols) - 1;       % 生成0/1序列
txSym = real(pammod(data, 2));           % 映射为±1（实部）

% 上采样并成形（无初始延迟）
txUpSequence = upsample(txSym, L, 0);    % 在符号间插入L-1个零

% 匹配滤波器输出（含成形滤波）
mfOutput = conv(r_p, txUpSequence);

% 导数匹配滤波器输出（仅MLTED使用）
dmfOutput = conv(conv(htx, txUpSequence), dmf);

%% 定时误差范围定义
tau = 0;  % 假设真实定时偏移为0（基准点）
tauEstVec = -(L/2):(L/2);      % 接收机估计的偏移范围：-L/2 ~ L/2 采样点
tauErrVec = tau - tauEstVec;   % 实际误差 = 真实偏移 - 估计偏移
normTauE = tauErrVec / L;      % 归一化为符号周期单位（-0.5 ~ 0.5）

%% S曲线计算（遍历所有定时误差）
g = zeros(1, length(tauErrVec));  % 初始化S曲线存储数组

for i = 1:length(tauEstVec)
    % 当前估计的定时偏移
    tauEst = tauEstVec(i);      % 接收机假设的定时偏移（采样点）
    tauErr = tau + tauEst;      % 实际定时误差

    % 理想符号采样点索引（补偿滤波器延迟）
    idealStrobeIdx = (rcDelay * L) + (1:L:(nSymbols * L));
    idealStrobeIdx_2sps = (rcDelay * L) + (1:L/2:(nSymbols * L));
    idealStrobeIdx_4sps = (rcDelay * L) + (1:L/4:(nSymbols * L));
    % 实际采样点索引（引入估计偏移）
    strobeIdx = idealStrobeIdx + tauEst;
    strobeIdx_2sps = idealStrobeIdx_2sps + tauEst;
    strobeIdx_4sps = idealStrobeIdx_4sps + tauEst;
    % 从匹配滤波器输出中提取符号值
    rxSym = mfOutput(strobeIdx);
    rxSym_2sps = mfOutput(strobeIdx_2sps);
    rxSym_4sps = mfOutput(strobeIdx_4sps);

    %% 根据TED类型计算误差信号
    switch (TED)
        case 'MLTED' % 最大似然TED（需导数滤波器）
            e = dmfOutput(strobeIdx) .* txSym;  % 导数采样 × 已知符号

        case 'ELTED' % 早-迟TED
            earlyIdx = strobeIdx + L/2;        % 超前L/2采样点
            lateIdx  = strobeIdx - L/2;        % 滞后L/2采样点
            e = txSym .* (mfOutput(earlyIdx) - mfOutput(lateIdx)); % 差值 × 符号

        case 'ZCTED' % 过零TED
            zcStrobeIdx = strobeIdx + L/2;     % 过零点采样位置
            zcSamples = mfOutput(zcStrobeIdx); % 提取过零点采样
            % 过零点采样 × 符号差分（前符号 - 当前符号）
            e = zcSamples(1:end-1) .* (txSym(1:end-1) - txSym(2:end));

        case 'GTED' % Gardner TED（非数据辅助）
            zcStrobeIdx = strobeIdx(2:end) - L/2;   % 符号间过零点位置  % promt 当前
            prevStrobeIdx = strobeIdx(1:end-1);     % 前符号位置        % early 较早
            currStrobeIdx = strobeIdx(2:end);       % 当前符号位置      % late  较晚
            % 过零点采样 × (前符号 - 当前符号)
            e = mfOutput(zcStrobeIdx) .* ...
                (mfOutput(prevStrobeIdx) - mfOutput(currStrobeIdx)); % 也是符合公式         
            %             e = real(conj(mfOutput(zcStrobeIdx)) .* ...
            %                 (mfOutput(currStrobeIdx) - mfOutput(prevStrobeIdx)));   % 这种TED的算法可能更为正确
        case 'MMTED' % Mueller-Muller TED
            prevStrobeIdx = strobeIdx(1:end-1);     % 前符号位置
            currStrobeIdx = strobeIdx(2:end);       % 当前符号位置
            % 前符号 × 当前采样 - 当前符号 × 前采样
            e = txSym(1:end-1) .* mfOutput(currStrobeIdx) - ...
                txSym(2:end)   .* mfOutput(prevStrobeIdx);

            %e_lee =  (mfOutput(prevStrobeIdx)+1j*mfOutput(currStrobeIdx))...
             %   .*( conj(mfOutput(prevStrobeIdx)) + 1j*conj(mfOutput(currStrobeIdx)) ) ;

        case 'GOTED'
            % 对信号进行分块处理
            NFFT = 1024;

            % --- 步骤1: 信号零填充 (Zero Padding) ---
            lpad = [rxSym_2sps, zeros(1, NFFT)];  % 问题出在这，应该是下采样后的信号进行

            % 如果mfOutput是列向量，则使用: lpad = [sigRxOffset; zeros(NFFT, 1)];

            % --- 步骤2: 计算块的数量 ---
            numBlocks = floor(length(lpad) / NFFT); % 使用floor确保是整数

%             NFFT=length(rxSym_2sps);
%             lpad = rxSym_2sps;
%             numBlocks=1;
            % --- 步骤3: 初始化误差存储 ---
            errors = zeros(numBlocks, 1); % 预分配列向量，提高效率

            % --- 步骤4: 循环处理每个块 ---
            for blockIdx = 1:numBlocks
                % --- 步骤4.1: 提取块 ---
                startIdx = (blockIdx - 1) * NFFT + 1;
                endIdx = startIdx + NFFT - 1;

                blk = lpad(startIdx:endIdx); % 提取当前块

                % --- 步骤4.2: 计算FFT ---
                blkFFT = fft(blk);

                % --- 步骤4.3: 计算Godard TED误差 ---
                error_val = sum(imag(blkFFT(1:NFFT/2).*conj(blkFFT(NFFT/2+1:end))));

                % 存储误差
                errors(blockIdx) = error_val/NFFT;
            end
            e=(errors);

           
        case 'GOTED_Mod'
            % 对信号进行分块处理
            NFFT = 2048;

            % --- 步骤1: 信号零填充 (Zero Padding) ---
            lpad = [rxSym_2sps, zeros(1, NFFT)];  % 问题出在这，应该是下采样后的信号进行

            % 如果mfOutput是列向量，则使用: lpad = [sigRxOffset; zeros(NFFT, 1)];

            % --- 步骤2: 计算块的数量 ---
            numBlocks = floor(length(lpad) / NFFT); % 使用floor确保是整数

            % --- 步骤3: 初始化误差存储 ---
            errors = zeros(numBlocks, 1); % 预分配列向量，提高效率

            % --- 步骤4: 循环处理每个块 ---
            for blockIdx = 1:numBlocks
                % --- 步骤4.1: 提取块 ---
                startIdx = (blockIdx - 1) * NFFT + 1;
                endIdx = startIdx + NFFT - 1;

                blk = lpad(startIdx:endIdx); % 提取当前块

                % --- 步骤4.2: 计算FFT ---
                blkFFT = fft(blk);

                % --- 步骤4.3: 计算Godard TED误差 ---
                %                 error_val = sum(imag(blkFFT(1:NFFT/2).*conj(blkFFT(NFFT/2+1:end))));
                error_val= (1/(2*pi)) * angle( sum( blkFFT(1:NFFT/2) .* conj( blkFFT(NFFT/2+1:end) ) ) );
                % 存储误差
                errors(blockIdx) = error_val/NFFT;
            end
            e=(errors);

        case 'Lee'
            % Lee 时域实现/ 频域实现
            %  zcStrobeIdx = strobeIdx(2:end) - L/2;
            %  prevStrobeIdx = strobeIdx(1:end-1);     % 当前符号位置
            %  currStrobeIdx = strobeIdx(2:end);       % 后一符号位置
            %  n=2:length(strobeIdx);
            %  phi_a=exp(-1j*n*pi);
            %  phi_b=exp(-1j*(n-1/2)*pi);
            %  误差
            %  e1 = sum(abs(currStrobeIdx).^2.*phi_a)+...
            %  sum(real(conj(prevStrobeIdx).*currStrobeIdx).*phi_b);
            %  e  =1/(2*pi)*(angle((e1)));
            % Lee 时域实现
            % 块的长度
            NFFT = length(rxSym_2sps);
            g_bias=1; % 偏执项
            % --- 步骤2: 计算块的数量 ---
            numBlocks = floor(length(rxSym_2sps) / NFFT); % 使用floor确保是整数
            % --- 步骤3: 初始化误差存储 ---
            errors = zeros(numBlocks, 1); % 预分配列向量，提高效率

            for blockIdx=1:numBlocks
                px=rxSym_2sps((blockIdx-1)*NFFT+1:blockIdx*NFFT);
                L1 = length(px);
                n = 1:L1;
                % cosine part
                ex1 = (-1).^(n-1);          % ex1 = exp(-1j.*(ii-1).*pi);
                sum_1 = sum( abs(px).^2 .* ex1 );

                % sine part
                ex2 = 1j * (-1).^(n-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
                xh = px(2:end);
                xx = px(1:end-1);
                ex2 = ex2(1:end-1);
                sum_2 = sum( real(conj(xx).*xh) .* ex2 );
                
                ex3 = exp(-1j.*(n-2).*pi);
                ex3 = ex3(1:end-1);
                sum_3 = ( sum ( (xx+1j*xh).* (conj(xx)+1j*conj(xh)) .* ex3) );

                % with biasing
                s = g_bias*sum_1 + sum_2;
                s = conj(s);

                temp = 1/2/pi * angle(s);
                temp_1 =1/2/pi * angle(sum_3);
                errors(blockIdx)=temp_1;
            end
            e=(errors);
        case 'SLN'
            % 块的长度
            NFFT = length(rxSym_4sps);
            % --- 步骤2: 计算块的数量 ---
            numBlocks = floor(length(rxSym_4sps) / NFFT); % 使用floor确保是整数
            % --- 步骤3: 初始化误差存储 ---
            errors = zeros(numBlocks, 1); % 预分配列向量，提高效率
            for blockIdx=1:numBlocks
                px=rxSym_4sps((blockIdx-1)*NFFT+1:blockIdx*NFFT);
                N = length(px);
                k = 1:N;
                ex = exp(-1j.*(k-1).*pi./2);
                s = sum( abs(px).^2 .* ex );
                temp = 1/2/pi*angle(s);
                errors(blockIdx)=temp;
            end
            e=errors;
    end

    %% 计算当前误差对应的S曲线值
    idx = tauErrVec == tauErr;       % 定位当前误差在向量中的索引
    g(idx) = mean(e);               % 误差信号取平均作为S曲线值

end

end % 函数结束
