function [ySym, H, errSq, HHist] = coreAdaptEq2(x, dx, SpS, H, nSymb, ...
                                                stepSize, lambdaRLS, ...
                                                eqLen, storeCoeff, alg, constSymb)
%COREADAPTEQ 自适应均衡核心处理函数
%
%   [ySym, H, errSq, HHist] = coreAdaptEq(x, dx, SpS, H, nSymb, ...
%                                         stepSize, lambdaRLS, ...
%                                         eqLen, storeCoeff, alg, constSymb)
% 输入
%   x         : L×N 复数接收信号（L 样点，N 模式）
%   dx        : nSymb×N 参考符号（训练序列）
%   SpS       : 每符号样点数
%   H         : N²×eqLen 当前抽头矩阵
%   nSymb     : 期望输出符号数
%   stepSize  : 算法步长 μ
%   lambdaRLS : RLS 遗忘因子（仅 RLS 类算法使用）
%   eqLen     : 均衡器长度（抽头数）
%   storeCoeff: 'true'|'false' 是否记录抽头历史
%   alg       : 算法标识 'nlms'|'cma'|'dd-lms'|'rde'|'da-rde'|'rls'|'dd-rls'|'static'
%   constSymb : 星座向量，部分算法需要
%
% 输出
%   ySym  : nSymb×N 均衡后符号
%   H     : N²×eqLen 更新后抽头矩阵
%   errSq : N×nSymb 每符号各模式瞬时误差平方
%   HHist : 抽头历史记录（storeCoeff=='true' 时有效）

%==========================================================================
%  0. 输入保护与常量
%==========================================================================
nModes = size(x, 2);
assert(size(dx,1)>=nSymb && size(dx,2)==nModes, '参考序列尺寸不符');
assert(size(H,1)==nModes^2 && size(H,2)==eqLen, '抽头矩阵尺寸不符');

%==========================================================================
%  1. 初始化
%==========================================================================
ySym      = zeros(nSymb, nModes);
errSq     = zeros(nModes, nSymb);
outEq     = zeros(nModes, 1);
idxTaps   = 1:eqLen;

% 抽头历史
if strcmp(storeCoeff, 'true')
    HHist = zeros(nModes^2, eqLen, nSymb);
else
    HHist = zeros(nModes^2, eqLen, 1);
end

% RLS 协方差初始化
if strcmp(alg, 'rls') || strcmp(alg, 'dd-rls')
    S = kron(eye(nModes), eye(eqLen));
end

% 算法相关常量
switch alg
    case {'cma', 'rde'}
        Rcma = complex(mean(abs(constSymb).^4) / mean(abs(constSymb).^2)* ones(1, nModes));
        Rrde = unique(abs(constSymb));
end

%==========================================================================
%  2. 逐符号处理
%==========================================================================
for sym = 1:nSymb
    % 2-1 计算当前输入窗
    idxSample = idxTaps + (sym-1)*SpS;      % 样点索引
    outEq(:)  = 0;
    
    % 2-2 MIMO 卷积
    for m = 1:nModes
        inSig = x(idxSample, m);            % 第 m 模式输入
        outEq = outEq + H((1:nModes) + (m-1)*nModes, :) * inSig;
    end
    ySym(sym, :) = outEq.';                 % 1×N 输出
    
    % 2-3 按算法更新抽头
    switch alg
        case 'nlms'
            [H, errSq(:,sym)] = nlmsUp(x(idxSample,:), dx(sym,:), outEq, stepSize, H, nModes);
        case 'cma'
            [H, errSq(:,sym)] = cmaUp(x(idxSample,:), Rcma, outEq, stepSize, H, nModes);
        case 'dd-lms'
            [H, errSq(:,sym)] = ddlmsUp(x(idxSample,:), constSymb, outEq, stepSize, H, nModes);
        case 'rde'
            [H, errSq(:,sym)] = rdeUp(x(idxSample,:), Rrde, outEq, stepSize, H, nModes);
        case 'da-rde'
            [H, errSq(:,sym)] = dardeUp(x(idxSample,:), dx(sym,:), outEq, stepSize, H, nModes);
        case 'rls'
            [H, S, errSq(:,sym)] = rlsUp(x(idxSample,:), dx(sym,:), outEq, lambdaRLS, H, S, nModes);
        case 'dd-rls'
            [H, S, errSq(:,sym)] = ddrlsUp(x(idxSample,:), constSymb, outEq, lambdaRLS, H, S, nModes);
        case 'static'
            errSq(:,sym) = 0;               % 无更新
        otherwise
            error('未知算法: %s', alg);
    end
    
    % 2-4 记录抽头
    if strcmp(storeCoeff, 'true')
        HHist(:,:,sym) = H;
    else
        HHist(:,:,1) = H;
    end
end
end % coreAdaptEq