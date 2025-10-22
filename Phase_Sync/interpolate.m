%% 插值函数
function [xI] = interpolate(method, x, m_k, mu, b_mtx, poly_f)
%根据样本向量x计算插值结果xI
%
% 参数说明:
%     method -> 插值方法: 多相滤波(0), 线性(1), 二次(2), 三次(3)
%     x      -> 用于计算插值的样本向量，包含基准点和周围的样本
%     m_k    -> 基准点索引，插值点之前的索引
%     mu     -> 基准点索引与期望插值时刻之间的估计得到的——分数间隔
%     b_mtx  -> 多项式插值系数矩阵，用于(二次)或(三次)插值
%     poly_f -> 多相滤波器组，当使用多相插值器(方法0)时处理输入样本

% 当mu超出标准范围[0,1)时调整基准点。这一步仅在支持奇数过采样率时需要，
% 此时会在原始mu估计上添加±0.5的偏移。相反，对于偶数过采样率，mu定义上就在[0,1)范围内
if (mu < 0)
    m_k = m_k - 1;   % mu小于0时，基准点减1
    mu = mu + 1;     % mu值加1，使其落入[0,1)
elseif (mu >= 1)
    m_k = m_k + 1;   % mu大于等于1时，基准点加1
    mu = mu - 1;     % mu值减1，使其落入[0,1)
end
assert(mu >= 0 && mu < 1);  % 确保mu在[0,1)范围内

switch (method)
    case 0 % 多相插值器
        % 使用mu选择多相子滤波器。使用floor操作确保得到的分支(polyBranch)
        % 始终在有效范围[1, polyInterpFactor]内，因为mu在[0,1)范围内。
        % 具体来说，添加第一个子滤波器的移位版本，即与第一个子滤波器相同但延迟缩短一个采样周期的子滤波器
        % 子滤波器相位已经非常接近。
        polyInterpFactor = size(poly_f, 1);  % 多相插值因子，即子滤波器数量
        polyBranch = floor(polyInterpFactor * mu) + 1;  % 计算选中的子滤波器分支
        polySubfilt = poly_f(polyBranch, :);  % 获取选中的子滤波器系数
        N = length(polySubfilt);              % 子滤波器长度
        % 用子滤波器与对应区间的样本卷积计算插值结果
        xI = polySubfilt * x((m_k - N + 1) : m_k);

    case 1 % 线性插值器
        % 线性插值公式：mu权重的下一个样本加上(1-mu)权重的当前样本
        xI = mu * x(m_k + 1) + (1 - mu) * x(m_k);

    case 2 % 二次插值器
        v_l = x(m_k - 1 : m_k + 2).' * b_mtx;  % 计算系数向量
        xI = (v_l(3) * mu + v_l(2)) * mu + v_l(1);  % 二次多项式计算

    case 3 % 三次插值器
        v_l = x(m_k - 1 : m_k + 2).' * b_mtx;  % 计算系数向量
        % 三次多项式计算
        xI = ((v_l(4) * mu + v_l(3)) * mu + v_l(2)) * mu + v_l(1);
end
end