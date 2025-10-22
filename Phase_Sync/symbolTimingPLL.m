function [xI]=symbolTimingPLL(TED, intpl, L, mfIn, mfOut, K1, K2, ...
    const, Ksym, rollOff, rcDelay, debug_s, debug_r,debug_r_fast,type)
% 时域方面  symbolPLL  算法
%  参数说明：
% TED：定时误差检测器类型: 'MLTED'/'ELTED'/'ZCTED'/'GTED'/'MMTED'/'SLN'/'Lee'
% intpl
% 误差检测器包括:
%     - Maximum-likelihood TED (MLTED);
%     - Early-late TED (ELTED);
%     - Zero-crossing TED (ZCTED);
%     - Gardner TED (GTED);
%     - Mueller-Muller TED (MMTED).
% 插值算法包括:
%     - Polyphase;
%     - Linear;
%     - Quadratic;
%     - Cubic.
% 输入信号
% mfin 匹配滤波器输入信号（多相插值时使用）
% mfout 匹配滤波器输出信号（非多相插值时使用）
% K1/K2   - PI控制器比例/积分增益
% const   - 星座图（用于符号判决）
% Ksym    - 符号归一化因子（接收信号需除以该因子后判决）
% rollOff - 升余弦滚降因子（0~1）
% rcDelay - 升余弦滤波器时延（符号数，通常为Tx/Rx滤波器总时延）
% debug_s - 静态调试绘图开关（0关闭，1开启）
% debug_r - 实时调试示波器开关（0关闭，1开启）
% 输出:
% xI      - 同步后的符号序列（列向量）


%% 数据预处理；参数设定
% 参数校验与默认值设置
if (nargin < 12)
    debug_s = 0; % 默认关闭静态调试绘图
end
if (nargin < 13)
    debug_r = 0; % 默认关闭实时调试示波器
end

% 调制阶数提取，输入信号格式化，输入信号选择
M = numel(const); % 根据星座图获取调制阶数


% 将输入信号转为列向量（确保后续处理一致性）
if (size(mfIn, 1) == 1)
    mfIn = mfIn(:);
end
if (size(mfOut, 1) == 1)
    mfOut = mfOut(:);
end

% 输入信号源（多相插值用MF输入，其他用MF输出）
if (intpl == 0)
    inVec = mfIn;   % 多相插值：直接处理未滤波信号（联合完成MF）
else
    inVec = mfOut;  % 其他插值：需外部预先完成MF
end


%% 实时调试工具配置（星座图+定时误差监视）
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

% Time scope used to debug the fractional error
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
    bufferLength = 1e4;
    waveformBuffer = zeros(1, bufferLength);
    currentIndex = 0;
    % 设置坐标轴范围
    xlim(hAx, [1 bufferLength]);  % 时间跨度1e4个点
    ylim(hAx, [-1 1]);   % Y轴范围[-1,1]

    % 初始化波形对象 - 使用初始点避免空对象
    hLine = plot(hAx, 1, 0, 'b');  % 先绘制第一个点

end

%% 快速调试图
if (debug_r_fast)
    hScope = comm.ConstellationDiagram(...
        'SymbolsToDisplaySource', 'Property',...
        'SamplesPerSymbol', 1, ...
        'MeasurementInterval', 256, ...
        'ReferenceConstellation', ...
        const);
    hScope.XLimits = [-1.5 1.5]*max(real(const));
    hScope.YLimits = [-1.5 1.5]*max(imag(const));
end


if (debug_r_fast)
    hTScopeCounter = dsp.TimeScope(...
        'Title', 'Fractional Inverval', ...
        'NumInputPorts', 1, ...
        'ShowGrid', 1, ...
        'ShowLegend', 1, ...
        'BufferLength', 1e5, ...
        'TimeSpanOverrunAction', 'Scroll', ...
        'TimeSpan', 1e4, ...
        'TimeUnits', 'None', ...
        'YLimits', [-1 1]);
end

%% 多相插值准备
% Polyphase filter bank
if (intpl == 0)
    polyInterpFactor = 128; % 多相插值因子（细分128个相位）
    % 设计高倍过采样升余弦滤波器（联合实现插值与匹配滤波）
    interpMf = sqrt(polyInterpFactor) * ...
        rcosdesign(rollOff, rcDelay, L * polyInterpFactor);
    % 多相分解：将长滤波器拆分为polyInterpFactor个子滤波器
    polyMf = polyDecomp(interpMf, polyInterpFactor);
    assert(size(polyMf, 1) == polyInterpFactor);

    % 生成多相导数滤波器组（用于MLTED）
    polyDMf = zeros(size(polyMf));
    for i = 1:polyInterpFactor
        polyDMf(i, :) = derivativeMf(polyMf(i, :), L);
    end

    % 翻转滤波器系数（便于卷积计算
    polyMf = fliplr(polyMf);
    polyDMf = fliplr(polyDMf);
else
    polyMf = []; % dummy variable
end


%% Farrow结构系数矩阵（二次/三次插值）
if (intpl == 2)
    % 二次插值系数（α=0.5优化参数）
    alpha = 0.5;
    b_mtx = flipud(fliplr(...
        [+alpha,      -alpha, 0; ...
        -alpha, (1 + alpha), 0; ...
        -alpha, (alpha - 1), 1; ...
        +alpha,      -alpha, 0]));

    %      b_mtx = [0 1 0 0;...
    %         -alpha alpha-1 alpha+1 -alpha;...
    %         alpha -alpha -alpha alpha]; % b(i,l)
elseif (intpl == 3)
    % Table 8.4.2
    % 三次插值系数（标准Farrow系数）
    b_mtx = flipud(fliplr(...
        [+1/6,    0, -1/6, 0; ...
        -1/2, +1/2,   +1, 0; ...
        +1/2,   -1, -1/2, 1; ...
        -1/6, +1/2, -1/3, 0]));

    %     b_mtx = [0 1 0 0;...
    %         -1/3 -1/2 1 -1/6;...
    %         1/2 -1 1/2 0;...
    %         -1/6 1/2 -1/2 1/6];
else
    b_mtx = []; % 线性/多相插值无需系数
end

%%  TimingRecovery feedback
% MLTED专用导数滤波器生成（非多相插值时）
if (intpl ~= 0 && strcmp(TED, 'MLTED'))
    % 设计根升余弦匹配滤波器
    mf = rcosdesign(rollOff, rcDelay, L);
    % 导数滤波器（中心差分法）
    dmf = derivativeMf(mf, L);
    % 对MF输入信号进行导数滤波（得到dMF输出）
    dMfOut = filter(dmf, 1, mfIn);
end


if strcmp(type,'feedback')

    % Recovery Loop param
    % Constants
    nSamples = length(inVec);
    nSymbols = ceil(nSamples / L);

    % Preallocate
    xI = zeros(nSymbols, 1); % 输出符号缓存
    mu = zeros(nSymbols, 1); % 估计得到的分数间隔μ缓存
    v  = zeros(nSamples, 1); % PI控制器输出缓存
    e  = zeros(nSamples, 1); % TED定时误差缓存
    
    % NCO debug
    debug_nco = 0; 
    if debug_nco
        out_nco = 0;
        out_nco1 = 0;
    end
    % 状态变量初始化
    k      = 0; % 当前符号索引
    strobe = 0; % 插值触发标志（1表示需要插值）
    cnt    = 1; % 模1计数器（初始值设为1以对齐首符号）
    vi     = 0; % PI积分器状态
    last_xI = inVec(1); % 上一个插值符号（初始化为第一个采样）

    % 确定输入数组边界
    if (strcmp(TED, 'ELTED'))
        n_end = nSamples - L;
    elseif (intpl > 1)
        n_end = nSamples - 1;
    else
        n_end = nSamples;
    end

    % 中点偏移计算（兼容奇偶过采样因子L）
    midpointOffset = ceil(L / 2);       % 符号间中点偏移量（如L=4时为2）
    muOffset = midpointOffset - L/2;   % μ偏移补偿（奇数L时为0.5）

    % 循环起始点（多相插值时需考虑滤波器长度）
    if (intpl == 0)
        poly_branch_len = size(polyMf, 2);% 多相子滤波器长度
        n_start = max(1, poly_branch_len - L); % 确保足够采样数
        if (strcmp(TED, 'ELTED') || strcmp(TED, 'ZCTED') || ...
                strcmp(TED, 'GTED'))
            n_start = n_start + ceil(L/2); % 补偿零交叉点偏移
        end
    else
        n_start = 1; % 其他插值方法从第一样本点开始
    end

    % 主循环
    for n = n_start:n_end
        if strobe == 1
            % 进行插值操作
            xI(k) = interpolate(intpl, inVec, m_k, mu(k), b_mtx, polyMf);

            % Timing Error Detector:
            a_hat_k = Ksym * slice(xI(k) / Ksym, M); % Data Symbol Estimate
            switch (TED)
                case 'MLTED' % Maximum Likelihood TED
                    % dMF interpolant
                    if (intpl == 0)
                        xdotI = interpolate(intpl, mfIn, m_k, mu(k), ...
                            b_mtx, polyDMf);
                    else
                        xdotI = interpolate(intpl, dMfOut, m_k, mu(k), b_mtx);
                    end
                    % Decision-directed version of Eq. (8.98), i.e., Eq. (8.27)
                    % adapted to complex symbols:
                    e(n) = real(a_hat_k) * real(xdotI) + ...
                        imag(a_hat_k) * imag(xdotI);
                case 'ELTED' % Early-late TED
                    % Early and late interpolants
                    early_idx = m_k + midpointOffset;
                    late_idx = m_k - midpointOffset;
                    early_mu = mu(k) - muOffset;
                    late_mu = mu(k) + muOffset;
                    x_early = interpolate(intpl, inVec, early_idx, ...
                        early_mu, b_mtx, polyMf);
                    x_late = interpolate(intpl, inVec, late_idx, ...
                        late_mu, b_mtx, polyMf);
                    % Decision-directed version of (8.99), i.e., (8.34)
                    % adapted to complex symbols:
                    e(n) = real(a_hat_k) * (real(x_early) - real(x_late)) + ...
                        imag(a_hat_k) * (imag(x_early) - imag(x_late));
                case 'ZCTED' % Zero-crossing TED
                    % Estimate of the previous data symbol
                    a_hat_prev = Ksym * slice(last_xI / Ksym, M);

                    % Zero-crossing interpolant
                    zc_idx = m_k - midpointOffset;
                    zc_mu = mu(k) + muOffset;
                    x_zc = interpolate(intpl, inVec, zc_idx, zc_mu, ...
                        b_mtx, polyMf);

                    % Decision-directed version of (8.100), i.e., (8.37)
                    % adapted to complex symbols:
                    e(n) = real(x_zc) * ...
                        (real(a_hat_prev) - real(a_hat_k)) + ...
                        imag(x_zc) * (imag(a_hat_prev) - imag(a_hat_k));
                case 'GTED' % Gardner TED
                    % Zero-crossing interpolant, same as used by the ZCTED
                    zc_idx = m_k - midpointOffset;
                    zc_mu = mu(k) + muOffset;
                    x_zc = interpolate(intpl, inVec, zc_idx, zc_mu, ...
                        b_mtx, polyMf);

                    % Equation (8.101):
                    e(n) = real(x_zc) * (real(last_xI) - real(xI(k))) ...
                        + imag(x_zc) * (imag(last_xI) - imag(xI(k)));
                case 'MMTED' % Mueller and Müller TED
                    % Estimate of the previous data symbol
                    a_hat_prev = Ksym * slice(last_xI / Ksym, M);

                    % Decision-directed version of (8.102), i.e., (8.49)
                    % adapted to complex symbols:
                    e(n) = real(a_hat_prev) * real(xI(k)) - ...
                        real(a_hat_k) * real(last_xI) + ...
                        imag(a_hat_prev) * imag(xI(k)) - ...
                        imag(a_hat_k) * imag(last_xI);
            end

            % % 更新历史插值符号,为下一次循环
            last_xI = xI(k);

            % ================ 实时调试 ================
            if (debug_r)

                %  更新星座图（固定长度数据点）
                hScatter.XData = [hScatter.XData(2:end), real(xI(k))];
                hScatter.YData = [hScatter.YData(2:end), imag(xI(k))];

                % 控制更新速度（可调整）
                drawnow limitrate;

                % 更新实时μ值曲线
                currentIndex = currentIndex + 1;
                if currentIndex > bufferLength
                    % 当缓冲区满时，移除最旧的数据，添加新数据
                    waveformBuffer = [waveformBuffer(2:end), mu(k)];
                    xlim(hAx, [currentIndex-bufferLength+1 currentIndex]);
                else
                    % 缓冲区未满时，直接添加新数据
                    waveformBuffer(currentIndex) = mu(k);
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
            end

            if (debug_r_fast)
                step(hScope, xI(k)) % 更新星座图
                step(hTScopeCounter, mu(k)); % 更新μ值曲线
            end
        else
            % 非插值时刻误差置零
            e(n) = 0;
        end
        % ================ loop filter ================

        % ================ PI控制器——loop filter更新 ================
        vp   = K1 * e(n);        % 比例项
        vi   = vi + (K2 * e(n)); % 积分项（累积）
        v(n) = vp + vi;          % 控制器输出(loop filter环路滤波输出)
        loopIntOut(n)  = vi; % 积分输出

        % ================ 模1计数器控制 ================
        W = 1/L + v(n);     % 动态步长（基础步长+控制量）  
        % 这里，1/L 是标准采样间隔的倒数，v(n) 是PI控制器的输出，用于动态调整间隔。W 本质上表示当前周期长度（sample）。
        
        % W 相当于 loopFilterOut ， cnt 相当于 out_nco
        % out_nco 最后输出 应该为 out_nco - = loopFilterOut ， 再做mod-1 运算
        % 限定 out_nco的范围 为 （0-1）
        
        strobe = cnt < W;   % 检查相位累加器，即发生下溢（触发插值）
        
        if (strobe)
            k = k + 1; % 更新符号索引
            m_k = n; % 记录当前基带点索引
            mu(k) = cnt / W; % 计算分数间隔μ ，表示在当前间隔内相对于理想采样点的相位偏移。这个值常用于插值器（如立方插值）来精确计算采样值。
        end

        % 更新计数器
        cnt = mod(cnt - W, 1);  % 更新相位累加器 cnt。
        % cnt - W 表示减去当前间隔 W，然后取模 1，确保 cnt 保持在 [0, 1) 范围内。
        % 模拟了一个数控振荡器（NCO）的行为，其中 cnt 代表归一化相位。取模操作处理相位环绕，确保循环连续。

%         if debug_nco
%             % 正常nco输出
%             out_nco=out_nco-v(n);
%             nco_value(n)=out_nco;
%             % 改进后nco输出
%             out_nco1=out_nco1-W;
%             nco_value1(n)=out_nco1;
%             % 改进后的loopfilter
%             loopfilter(n) =  W;
%             
%             if out_nco1 > 1
%                 out_nco1 =out_nco1- 1;
%             elseif out_nco1 < -1
%                 out_nco1 = out_nco1+ 1;
%             end
%             nco_value_test(n)=out_nco1;
% 
%             % nco 取模运算
%             nco_mod(n)=mod(out_nco,1);
%             nco_mod1(n) = cnt;
%         end

    end % for start

%     if debug_nco
% 
%         figure;hold on;
%         plot(v);
%         plot(loopfilter);
%         legend('loop_filter','loop_filter_modfity')
% 
%         figure; hold on;
%         plot(nco_value);
%         plot(nco_value1);
%         legend('nco_out','nco_out_modfity')
% 
%         figure; hold on;
%         plot(nco_mod)
%         plot(nco_mod1)
%         legend('nco_out_mod','nco_out_mod_modfity')
% 
%         figure;
%         plot(nco_value_test)
%         legend('测试nco输出，不是用取余的方式')
%     end

    % 裁剪输出序列（去除未完成的插值）
    if (strobe) % ended on a strobe (before filling the k-th interpolant)
        xI = xI(1:k-1);
    else
        xI = xI(1:k);
    end



    %%  TimingRecovery feedford
elseif strcmp(type,'feedford')

    % 分组处理
    szBlock = 1024;
    % 在后面填充 szBlock 个零
    signalPad = [ inVec; zeros(szBlock,1)];

    % --- 步骤2: 计算块的数量 ---
    numBlocks = ceil(length(inVec) / szBlock); % 使用floor确保是整数
    % 计算输出信号长度
    output_length = numBlocks * (szBlock);
    xI = zeros(output_length, 1); % 预分配输出
    nco_values = zeros(length(signalPad), 1); % 存储NCO输出相位


    % --- 初始化DPLL状态 ---
    out_nco = 0;    % NCO当前输出相位
    intPart = 0;    % PI滤波器积分部分
    pointer = 1 ;   % 记录当前基带点索引
    %indB = 1 ;  % 索引值

    ind = 0:szBlock-1; % 区间索引
    %while pointer+szBlock+L-2 <= length(signalPad)
    for indB = 1:numBlocks

        % --- 1. 提取信号块 ---
        step_start = (indB - 1) * (szBlock) + 1; % MATLAB索引从1开始
        end_step = step_start + szBlock - 1;
        signal = signalPad(step_start:end_step);

        % ---TED ---
        switch (TED)
            case 'LeeTED'
                L1 = length(signal);
                n1 = 1:L1;
                % cosine part
                ex1 = (-1).^(n1-1).';          % ex1 = exp(-1j.*(ii-1).*pi);
                sum_1 = sum( abs(signal).^2 .* ex1);
                % sine part
                ex2 = 1j * (-1).^(n1-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
                xh = signal(2:end);
                xx = signal(1:end-1);
                ex2 = ex2(1:end-1).';
                sum_2 = sum( real(conj(xx).*xh) .* ex2 );
                % with bias (bias  =  1)
                s = sum_1 + sum_2;
                %s = conj(s);

                temp = 1/2/pi*angle(s);
                ted(indB)=temp;
            case 'SLNTED'
                % 需要四倍上采样
                assert(L>2);
                L1 = length(signal);
                n1 = 1:L1;
                ex = exp(-1j.*(n1-1).*pi./2).';
                s = sum( abs(signal).^2 .* ex );
                temp = 1/2/pi*angle(s);
                ted(indB)=temp;
        end

        % --- 5. PI 环路滤波器 ---
        intPart = intPart + K2 * ted(indB); % 积分部分更新
        propPart = K1 * ted(indB);          % 比例部分
        loopFilterOut(indB) = propPart + intPart;
        

        % =========================================================
        %  块间定时漂移修正：先算差值，再决定“跳点/等点”
        % =========================================================

%         if indB == 1
%             tau1 = 0;                    % 第一块没有“上一块”，设 0
%         else
%             tau1 = loopFilterOut(indB-1);           % 上一块的定时偏移
%         end
%         tau2 = loopFilterOut(indB);                 % 当前块的定时偏移
% 
%         % 把 tau 差值转成“整点跳/等” + 小数余量 nu        
%         % Nco control
%         if floor((0.5-tau2)*L)
%             if tau1-tau2>0.75
%                 pointer = pointer+1;    % skip one sample
%                 tau2 = tau2+(L-1)/L;
%                 out_nco = tau2*L;
%             elseif tau1-tau2<=0.25
%                 pointer = pointer-1;     % wait for one sample
%                 tau2 = tau2+1/L;
%                 out_nco = tau2*L;
%             elseif tau1-tau2<=0.75 && tau1-tau2>0.5
%                 pointer = pointer+2;     % skip 2 sample
%                 tau2 = tau2+(L-2)/L;
%                 out_nco = tau2*L;
%             elseif tau1-tau2<=0.5 && tau1-tau2>0.25
%                 pointer = pointer-2;     % wait for 2 sample
%                 tau2 = tau2+2/L;
%                 out_nco = tau2*L;
%             end
%         else
%             out_nco = tau2*L;
%         end
%         out_nco = mod(out_nco,1);
% 
%         % 更新PI环路滤波输出
%         loopFilterOut(indB)=tau2;


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


        for m_k = 1 : szBlock
            % ================ interplate ================
            % 可能存在问题， 第一个插值，没有足够的符号(待解决)【高维度插值的问题】
            % 解决方法1 ： 使用线性插值
            y_out(m_k+(indB-1)*szBlock,1) = interpolate(1, signalPad, m_k+(indB-1)*szBlock, out_nco, b_mtx, polyMf);
            %(indB-1)*szBlock
        end
        % 更新采样基点
        pointer = pointer+szBlock;
        % pl(indB) = pointer;
        % indB = indB+1;
        % --- 10. 记录NCO值 ---
        nco_values(step_start:step_start+szBlock-1) = out_nco;
        
    end
    % 输出向量
    % y_out=y_out(:);
    xI  =   y_out(1:length(inVec));
     


    % Recovery Loop param
%     % Constants
%     nSamples = length(inVec);
%     nSymbols = ceil(nSamples / L);
% 
%     % Preallocate
%     xI = zeros(nSymbols, 1); % 输出符号缓存
%     mu = zeros(nSymbols, 1); % 估计得到的分数间隔μ缓存
%     v  = zeros(nSamples, 1); % PI控制器输出缓存
%     e  = zeros(nSamples, 1); % TED定时误差缓存
% 
%     % 状态变量初始化
%     k      = 0; % 当前符号索引
%     strobe = 0; % 插值触发标志（1表示需要插值）
%     cnt    = 1; % 模1计数器（初始值设为1以对齐首符号）
%     vi     = 0; % PI积分器状态
%     last_xI = inVec(1); % 上一个插值符号（初始化为第一个采样）
%     % 确定数据起始位置
%     n_start = 1 ;
%     % 结束位置
%     if (intpl > 1)
%         n_end = nSamples - 1;
%     else
%         n_end = nSamples;
%     end % for intpl
% 
%     % start Timign
%     for n = n_start:n_end
%         % 信号分组(有待讨论)
%         if strobe == 1
%             xI(k)=inVec(m_k) ;
%             m_k_inter = m_k;
%             % Timing Error Detector:
%             switch (TED)
%                 case 'LeeTED'
%                     %                     L1 = length(inVec);
%                     %                     n1 = 1:L1;
%                     % cosine part
%                     %                     ex1 = (-1).^(n1-1).';          % ex1 = exp(-1j.*(ii-1).*pi);
%                     %                     sum_1 = sum( abs(inVec).^2 .* ex1);
%                     ex1 = (-1).^k.';
%                     sum_1 = abs(xI(k)).^2.*ex1 ;
%                     % sine part
%                     %                     ex2 = 1j * (-1).^(n1-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
%                     %                     xh = inVec(2:end);
%                     %                     xx = inVec(1:end-1);
%                     %                     ex2 = ex2(1:end-1).';
%                     %                     sum_2 = sum( real(conj(xx).*xh) .* ex2 );
%                     ex2 = 1j * (-1).^(k);
%                     sum_2 = real(conj(xI(k)).*last_xI) .* ex2;
%                     % with bias (bias  =  1)
%                     s = sum_1 + sum_2;
%                     s = conj(s);
% 
%                     temp = -1/2/pi*angle(s);
%                     e(n)=temp;
%                 case 'SLNTED'
%                     % 需要四倍上采样
%                     assert(length(L)>2);
%                     %                     N1 = length(inVec);
%                     %                     n1 = 1:N1;
%                     %                     ex = exp(-1j.*(n1-1).*pi./2);
%                     %                     s = sum( abs(inVec).^2 .* ex );
%                     ex = 1j .^ (-(k-1));          % 结果：1, -j, -1, j, 1, -j, …
%                     s = abs( xI(k).^2 .* ex );
%                     temp = 1/2/pi*angle(s);
%                     e(n)=temp;
%             end
% 
% 
%         else
%             e(n)=0;
%         end % for strob
%         % ================ PI控制器更新 ================
%         vp   = K1 * e(n);        % 比例项
%         vi   = vi + (K2 * e(n)); % 积分项（累积）
%         v(n) = vp + vi;          % 控制器输出
%         % ================ 模1计数器控制 ================
%         W = 1/L + v(n);     % 动态步长（基础步长+控制量）
% 
%         strobe = cnt < W;   % 检查是否下溢（触发插值）
%         pp(n)= double(strobe);
%         if (strobe)
%             k = k + 1; % 更新符号索引
%             m_k = n; % 记录当前基带点索引
%             mu(k) = cnt / W; % 计算分数间隔μ
%             if k ==1
%                 xI(k) =  inVec(m_k);
%             else
%                 % ================ interplate ================
%                 xI(k-1) = interpolate(intpl, inVec, m_k_inter, mu(k), b_mtx, polyMf);
%                 % % 更新历史插值符号,为下一次循环
%                 last_xI = xI(k-1);
%             end
% 
%         end
% 
%         % 更新计数器
%         cnt = mod(cnt - W, 1);
% 
%     end % for nstart
% 
% 
%     % 裁剪输出序列（去除未完成的插值）
%     if (strobe) % ended on a strobe (before filling the k-th interpolant)
%         xI = xI(1:k-1);
%     else
%         xI = xI(1:k);
%     end

end % for ford or back


%% 算法运行图例及参数变化图
if (debug_s)
    figure
    plot(e)
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
    title('TED output')

    figure
    plot(v)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(loopIntOut)
    title('PI-Popter Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    
    figure
    plot(mu, '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
end


end