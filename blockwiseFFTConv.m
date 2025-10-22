function y = blockwiseFFTConv(x, h, NFFT, freqDomainFilter)
%BLOCKWISEFFTCONV 重叠-保留频域块卷积（Overlap-Save）
%
%   y = blockwiseFFTConv(x, h, NFFT, freqDomainFilter)
%
% 输入
%   x                : L×1 实/复向量，待滤波信号
%   h                : M×1 实/复向量，FIR 冲激响应；或 NFFT×1 频响（见下方）
%   NFFT             : 可选，FFT 长度（>M）。[] 或缺失时自动取 2 的幂
%   freqDomainFilter : 可选，true 表示 h 已是频响；false（默认）表示时域冲激响应
%
% 输出
%   y : L×1 向量，滤波结果，长度与 x 严格一致，群延迟已切除

    %% 0. 参数解析与基本检查
    if nargin < 3 || isempty(NFFT)
        NFFT = 2^nextpow2(numel(h));    % 最小 2 的幂 ≥ M
    end
    if nargin < 4
        freqDomainFilter = false;
    end

    x = x(:);                           % 强制列向量
    h = h(:);
    Lx  = numel(x);
    M   = numel(h);
    D   = floor((M-1)/2);               % 群延迟
    if NFFT < M
        error('NFFT 必须 ≥ 滤波器长度 M');
    end
    L   = NFFT - M + 1;                 % 每块有效输出长度

    %% 1. 滤波器预处理：得到 NFFT 点频响 H
    if freqDomainFilter
        % h 已是频响 → 先 IFFT 回时域，再补零，再 FFT
        hTime = ifft(fftshift(h));
        hTime = [hTime; zeros(NFFT-M,1)];
        H = fft(hTime);
    else
        % h 是冲激响应 → 直接补零后 FFT
        hPad = [h; zeros(NFFT-M,1)];
        H = fft(hPad);
    end

    %% 2. 信号补零：尾部凑整 + 群延迟补偿
    numBlocks = ceil(Lx / L);
    padLen    = numBlocks*L - Lx;       % 凑整块
    xPad      = [x; zeros(padLen + D,1)];% 多补 D 点，后期一次性切掉
    xPad      = [zeros(M-1,1); xPad];   % 头部再补 M-1 个 0，供第一块重叠

    %% 3. 输出缓存
    yFull = zeros(size(xPad));          % 先存全部，再切

    %% 4. 重叠-保留主循环
    startIdx = 1;
    for blk = 1:numBlocks
        segment = xPad(startIdx+(0:NFFT-1));    % 取 NFFT 点（含旧样点）
        X = fft(segment);
        yBlk = ifft(X .* H);                    % 频域乘
        % 只留后半段 L 点
        yFull((blk-1)*L + (1:L)) = yBlk(M:NFFT);
        startIdx = startIdx + L;                % 滑窗 L 点
    end

    %% 5. 切除群延迟 & 尾部补零段，返回与原始 x 同长度
    y = yFull( (M-1+1+D) : (M-1+Lx+D) );
end