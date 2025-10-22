% ============== 偏振态跟踪算法实现 ============== %
clc;clear;close all;

%% 信号生成
coherentGenerationClock;
% 添加路径
addpath('Plot\')
addpath('Dsp\')
addpath('Phase_Sync\')
addpath('Fncs\')
%% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% 发射信号 
[txSig,qamSig]  =  Tx.dataOutput();

% 参考信号
label=qamSig;
% 数据性能起始位置
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）

%% Sop tracking
outSignal=Rx.addRsop(txSig);

outSignal_X = downsample(outSignal(:,1),Tx.TxPHY.sps);
outSignal_Y = downsample(outSignal(:,2),Tx.TxPHY.sps);
scatterplot(outSignal_Y)


if 0
% SOP tracking
switch lower(SOPMethod)
    case 'nebojsa'
        % 1. 设置处理块大小，最大不超过1024个样本
        BlkSz = min(2^10,size(x,1));
        % 2. 初始化偏振态参数向量[SOP_theta, SOP_phi]，通常代表琼斯矩阵中的角度
        SOP = [0 0];
        % 3. 设置算法步长（学习率），将0.2度转换为弧度
        mu = 0.2/180*pi;
        % 4. 核心迭代循环：进行1000次迭代以收敛到最佳偏振态估计
        for iter = 1:1000
            % 调用Nebojsa算法函数，输入当前数据块、块大小、上一次的SOP估计和步长
            % 输出：~（忽略校正后的信号），更新的SOP参数，代价函数值Kd
            [~, SOP(iter+1,:), Kd(iter)] = Nebojsa(x(1:BlkSz,:),BlkSz,SOP(iter,:),mu);
        end
        % 5. 迭代收敛后，使用最终（第1001个）的SOP参数来构建补偿用的琼斯矩阵 (J)
        J = [cos(SOP(end,1))*exp(-1i*SOP(end,2)/2) -sin(SOP(end,1))*exp(1i*SOP(end,2)/2);...
            sin(SOP(end,1))*exp(-1i*SOP(end,2)/2) +cos(SOP(end,1))*exp(1i*SOP(end,2)/2)];
        % 6. 应用琼斯矩阵对整个输入信号x进行偏振补偿
        %    先将x转置，与J相乘，再转置回来，得到补偿后的信号x
        x = (J*x.').';
        % 7. （注释掉的代码）用于绘制SOP参数的收敛轨迹，用于调试和可视化
        %         plot(SOP(:,1),SOP(:,2));
    otherwise
        % 如果SOPMethod不是'nebojsa'或未指定，则不进行偏振态跟踪与补偿
end


% 估计DGD
% DGD estimate
% 1. 为DGD估计设置一个新的块大小，通常比SOP跟踪的块更大（最大4096样本），以获得更稳定的统计结果
BlkSz = min(2^12,size(x,1));
% 2. 计算DGD估计值
DGDest = -1/2/pi*angle(ted_GodardFF(x(1:BlkSz,1)))+1/2/pi*angle(ted_GodardFF(x(1:BlkSz,2)));

end 


%%
% Nebojsa偏振态跟踪算法
function [y, paraout, Kd] = Nebojsa(x,bs,parain,mu)
X = fft(x(:,1))/sqrt(size(x,1));
Y = fft(x(:,2))/sqrt(size(x,1));
stepsz = [1 1;1 -1;-1 1;-1 -1]*mu; % 步长矩阵
para = [parain(1)+stepsz(:,1) parain(2)+stepsz(:,2)]; % 参数候选
for n = 1:size(para,1)
    % 应用琼斯矩阵变换
    Z1 = X*cos(para(n,1))*exp(-1i*para(n,2)/2)-Y*sin(para(n,1))*exp(1i*para(n,2)/2);
    Z2 = Z1;
    % 计算代价函数
    KdCdd(n) = real(Z1(1:bs/2).'*Z2(bs/2+1:bs)'.');
end
[Kd,ind] = max(KdCdd); % 选择最佳参数
paraout = para(ind,:);

% 构建最佳琼斯矩阵并应用
J = [cos(paraout(1))*exp(-1i*paraout(2)/2) -sin(paraout(1))*exp(1i*paraout(2)/2);...
    sin(paraout(1))*exp(-1i*paraout(2)/2) +cos(paraout(1))*exp(1i*paraout(2)/2)];
y = (J*x.').'; % 偏振补偿
end

% Lingchen偏振态跟踪算法
function [y, paraout, Kd] = Lingchen(x,bs,parain,mu)
X = fft(x(:,1))/sqrt(size(x,1));
Y = fft(x(:,2))/sqrt(size(x,1));
stepsz = [1 0;-1 0]*mu; % 步长矩阵（仅调整一个参数）
para = [parain(1)+stepsz(:,1) parain(2)+stepsz(:,2)]; % 参数候选
for n = 1:size(para,1)
    % 应用琼斯矩阵变换
    Z1 = +X*cos(para(n,1))*exp(+1i*para(n,2)/2)+Y*sin(para(n,1))*exp(+1i*para(n,2)/2);
    Z2 = -X*sin(para(n,1))*exp(-1i*para(n,2)/2)+Y*cos(para(n,1))*exp(-1i*para(n,2)/2);
    % 计算代价函数
    KdCdd(n) = abs(Z1(1:bs/2).'*Z1(bs/2+1:bs)'.')+abs(Z2(1:bs/2).'*Z2(bs/2+1:bs)'.');
end
[Kd,ind] = max(KdCdd); % 选择最佳参数
paraout = para(ind,:);

% 构建最佳琼斯矩阵并应用
J = [+cos(paraout(1))*exp(+1i*paraout(2)/2) +sin(paraout(1))*exp(+1i*paraout(2)/2);...
    -sin(paraout(1))*exp(-1i*paraout(2)/2) +cos(paraout(1))*exp(-1i*paraout(2)/2)];
y = (J*x.').'; % 偏振补偿
end