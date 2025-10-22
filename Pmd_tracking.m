% ============== PMD-模型实现 ============== %
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
%% PMD model —— 张晓光

% tau_mean = 3*1e-12; % average tau 3ps
% sigma = 0.01; % variance
% w=2*pi*193.5*1e12;
% % 输入信号为2*n形式
% outSignal=Rx.addPMD(txSig.',tau_mean,sigma,w);
% outSignal=outSignal.';
% outSignal_X = downsample(outSignal(:,1),Tx.TxPHY.sps);
% outSignal_Y = downsample(outSignal(:,2),Tx.TxPHY.sps);
% scatterplot(outSignal_Y)
%% PMD 简便实现模型
[f_pmd, t_pmd]=freq_time_set(size(txSig,1),fs);
w_pmd = 2*pi * f_pmd.';
% PDL损耗
meanPDLdB=1;             %引入PDL相关损耗
gama=10^(-meanPDLdB/20);    %注意在这里加入了DGD损耗
% DGD设置—— 以一个符号周期为标尺尺度
DGD_symbl = 1/fb;
meanDGD = DGD_symbl*0.4; % 平均的DGD设置为0.4个符号周期
PP1= exp(1j*w_pmd*meanDGD/2);  % 快轴相移，对角相，引入DGD
PP2 = 0;
% RSOP旋转
rsop_theta = pi/8;

% 将X和Y偏振信号转换到频域
sig_fft_x = fft(sigin(:,1));
sig_fft_y = fft(sigin(:,2));
% 将零频率分量移到频谱中心
sig_fft_x = fftshift(sig_fft_x);
sig_fft_y = fftshift(sig_fft_y);

% PMD中DGD应用
uux = gama.*PP1.* sig_fft_x + PP2.* sig_fft_y;
uuy = -conj(PP2).* sig_fft_x + conj(PP1).*sig_fft_y;
% 转换为时域
sig_x=ifft(fftshift(uux));            %%这里对应下文的re_ex
sig_y=ifft(fftshift(uuy));
% RSOP 旋转
mat_rot = [cos(rsop_theta),-sin(rsop_theta);...
           sin(rsop_theta),cos(rsop_theta)];
% SOP 应用
sig_x_out = mat_rot(1, 1) .* sig_fft_x + mat_rot(1, 2).* sig_fft_y;
sig_y_out = mat_rot(2, 1).* sig_fft_x + mat_rot(2, 2).* sig_fft_y;
% 返回处理后的信号
sigout(:,1) = sig_x_out;
sigout(:,2) = sig_y_out;

%% PMD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    pmd_config:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         dgd_manual: 1            # 是否设定一个跨段总差分群时延大小（1 - 是，0 - 否）
%         dgd_total: 0.2            # 单跨段差分群时延（ps），dgd_manual=1需设置此参数
%         psp_manual: 1             # 是否固定两个trunk之间主偏振态的旋转角度（1 - 是，0 - 否）,为1则角度为固定值pi/4
%  0 (默认)：随机模式。PSP 的方向在庞加莱球上均匀随机分布，模拟真实光纤中双折射轴的随机变化。
% 1：固定模式。使用固定的角度 phi 作为所有段的 PSP 方向。适用于简化测试或研究固定双折射的影响
%         pmd_coeff_random: 0        # 每个trunk的pmd系数是否随机（1 - 是，0 - 否）【pmd_coeff固定 或者从麦克斯韦分布中采样】
%         pmd_dz_random: 0           # 各trunk长度是否随机（1 - 是，0 - 否），为0时各trunk等长
%         pmd_coeff: 0.05            # 单个跨段的pmd系数（ps.km^(-1/2)），dgd_manual=0需设置此参数，默认值0.05
%         pmd_trunk_num: 80          # 单个跨段的trunk数目，默认值80

% dgd_rms = sqrt(3 * pi / 8 * pmd_coeff.^2)
% if pmd_coeff_random:
%     # sample the pmd coefficient from Maxwellian distribution
%     vx = np.random.normal(loc = 0, scale = sqrt(dgd_rms**2/3))
%     vy = np.random.normal(loc = 0, scale = sqrt(dgd_rms**2/3))
%     vz = np.random.normal(loc = 0, scale = sqrt(dgd_rms**2/3))
%     pmd_coeff = np.sqrt(vx**2 + vy**2 + vz**2)
% else:
%     pmd_coeff = pmd_coeff (固定)

% 相关长度:
% l_corr = span_len/pmd_trunk_num (每一段分为多少来模拟PMD)

% SOP的旋转角度
%         if psp_manual:
%             固定角度
%             self.psp_theta   = phi * ones(pmd_trunk_num)
%             self.psp_phi     = phi * ones(pmd_trunk_num)
%         else:
%             self.psp_theta = random.rand(pmd_trunk_num) * 2 * np.pi - np.pi           # 均匀分布[0-1） azimuth: uniform R.V.
%             self.psp_phi   = 0.5 * arcsin(np.random.rand(pmd_trunk_num) * 2 - 1)
%                 # uniform R.V. over the Poincare sphere


% GVD参数：
%  计算角频率数组 (Angular Frequency Array)
%  2 * pi * f: 将频率f (Hz) 转换为角频率 ω (rad/s)
[f, t]=freq_time_set(size(txSig,1),fs);
w = 2*pi * f.';
%  计算色散算子在频域的相位因子 (Phase Factor for Chromatic Dispersion in Frequency Domain)
%  该因子是光脉冲在光纤中传播常数β(ω)的泰勒展开式，用于描述色散（CD）效应
%  在频域中，色散效应等价于乘以一个相位因子 exp(-j * phase_factor_freq * dz)
% phase_factor_freq = beta0+ beta1 * (w) + beta2 * (w.^ 2) / 2+ beta3 * (w.^ 3) / 6 ;
phase_factor_freq=0;
%  零阶项: β0, 代表相位常数，通常合并到载波中，这里为完整性而保留
%  一阶项: β1 * ω, 代表群时延 (Group Delay)，决定脉冲的传输速度
%  二阶项: (β2 * ω^2) / 2, 代表群速度色散 (GVD - Group Velocity Dispersion)，是导致脉冲展宽的主因
%  三阶项: (β3 * ω^3) / 6, 代表三阶色散 (Third-Order Dispersion) 或色散斜率 (Dispersion Slope)

dz=0.5; %0.5km
pmd_trunk_num=80;   %单个跨段的trunk数目，默认值80
l_corr=dz/pmd_trunk_num; % 相关长度
% 旋转角度设置
psp_manual = 1;
if psp_manual
    psp_theta = pi/4 * ones(pmd_trunk_num, 1);
    psp_phi = pi/4 * ones(pmd_trunk_num, 1);
else
    % 均匀分布在[-π, π)区间
    psp_theta = rand(pmd_trunk_num, 1) * 2 * pi - pi;
    % 在Poincare球面上的均匀分布
    psp_phi = 0.5 * asin(rand(pmd_trunk_num, 1) * 2 - 1);
end

% 参数设置：
kwargs=struct();
kwargs.trunk_list= l_corr:l_corr:dz;  % PMD小段列表  设置为相同长度，或者随机长度，长度数量确定 【可进阶】
kwargs.psp_theta=psp_theta ;  % 主态方向角度θ
kwargs.psp_phi=psp_phi ;  % 主态方向角度φ

pmd  = 1 ; % pmd仿真开关
pmd_coeff = 0.05; %单个跨段的pmd系数（ps.km^(-1/2)），dgd_manual=0需设置此参数，默认值0.05


if 0
    dgd_manual=1 ;
    dgd_total = 0.2 ;
    pmd_coeff_random=1;
    if dgd_manual
        pmd_coeff= dgd_total/dz;
    else
        dgd_rms = sqrt(3 * pi / 8 * pmd_coeff.^2);
        if pmd_coeff_random
            vx = normrnd(0, sqrt(dgd_rms.^2 / 3));
            vy = normrnd(0, sqrt(dgd_rms.^2 / 3));
            vz = normrnd(0, sqrt(dgd_rms.^2 / 3));
            pmd_coeff = sqrt(vx^2 + vy^2 + vz^2);
        else
            pmd_coeff =pmd_coeff ;
        end % for pmd_coeff_random
    end % for dgd_manual
end


% 每个分段的PMD参数： pmd_per_trunk
pmd_per_trunk =  pmd_coeff / sqrt(pmd_trunk_num); % PMD coefficient per trunk.

% 计算差分相移的因子，与（色散）对应，决定DGD的大小
differential_factor_freq = 0.5 * w * pmd_per_trunk;

% PMD建模
sigout_pdm = Pmd_model(txSig, dz, phase_factor_freq, pmd, differential_factor_freq, kwargs);



function sigout = Pmd_model(sigin, dz, phase_factor_freq, pmd, differential_factor_freq, kwargs)
% 计算双偏振信号的色散效应
%
% 此函数计算双偏振信号的色散效应，包括群速度色散(GVD)和偏振模色散(PMD)，其中PMD是可选的。
% 首先对色散算子和信号进行傅里叶变换，然后在频域中使用矩阵乘法添加色散效应。
%
% 输入参数:
%   sigin : 元胞数组
%       传输信号，包含X和Y两个偏振分量。
%   dz : 浮点数
%       本次计算需要处理的色散长度。
%   phase_factor_freq : 数组
%       色散算子的傅里叶变换，用于计算GVD效应。
%   pmd : 整数, {0,1}
%       确定仿真中是否考虑PMD效应。
%   differential_factor_freq : 数组, 可选
%       差分群延时(DGD)算子的傅里叶变换，用于计算PMD效应。
%   data_mode : 字符串, {'numpy','tensor'}, 可选
%       保留参数，为了与Python版本兼容，在MATLAB中通常使用默认处理方式。
%   kwargs : 结构体
%       与PMD效应仿真相关的参数。
%
% 输出参数:
%   sigout : 元胞数组


% 将X和Y偏振信号转换到频域
sig_fft_x = fft(sigin(:,1));
sig_fft_y = fft(sigin(:,2));
% 将零频率分量移到频谱中心
sig_fft_x = fftshift(sig_fft_x);
sig_fft_y = fftshift(sig_fft_y);

if pmd  % 如果启用PMD仿真
    % 从结构体中获取PMD相关参数
    %     trunk_idx_list = kwargs.trunk_idx_list;  % PMD小段索引列表   长度数量确定
    trunk_list = kwargs.trunk_list;  % PMD小段列表  设置为相同长度，或者随机长度，长度数量确定 【可进阶】
    psp_theta = kwargs.psp_theta;  % 主态方向角度θ
    psp_phi = kwargs.psp_phi;  % 主态方向角度φ

    % 定义Pauli矩阵
    sig0 = eye(2);  % 单位矩阵
    sig2 = [0, 1; 1, 0];  % Pauli σ2矩阵
    sig3i = [0, 1; -1, 0];  % Pauli σ3矩阵的虚数倍

    prop_dz = 0;  % 已处理的累计长度

    % 按顺序将信号通过PMD小段
    for i = 1:length(trunk_list)
        pmd_dz = trunk_list(i);  % 当前PMD小段的长度

        % 构建旋转矩阵的一部分（θ旋转）
        mat_theta = cos(psp_theta(i)) * sig0 ...
            - sin(psp_theta(i)) * sig3i ;  % 正交矩阵  % 不确定是否需要写成复数形式

        % 构建旋转矩阵的另一部分（φ旋转）
        mat_epsilon = cos(psp_phi(i)) * sig0 ...
            + 1j * sin(psp_phi(i)) * sig2;  % 正交矩阵

        % 组合成完整的旋转矩阵：用于将信号变换到主态(PSP)基
        mat_rot = mat_theta * mat_epsilon;  % 得到RSOP矩阵
        mat_rot_conj = conj(mat_rot);  % 旋转矩阵的共轭

        % 注意：设A=[Ax;Ay]为电场，则matR*D*matR'*A是线性PMD步骤，其中D是DGD矩阵

        % 1> 通过乘以旋转矩阵R的共轭，将信号变换到PSP基上
        uux = mat_rot_conj(1, 1) * sig_fft_x + mat_rot_conj(2, 1) * sig_fft_y;
        uuy = mat_rot_conj(1, 2) * sig_fft_x + mat_rot_conj(2, 2) * sig_fft_y;

        % 2> 应用双折射、DGD和GVD：全部包含在对角矩阵D中
        gvd_beta = phase_factor_freq * pmd_dz;  % GVD β因子
        pmd_beta = sign(pmd_dz) * differential_factor_freq * sqrt(abs(pmd_dz));  % 差分β因子
        % 注意：dzb(k)/brf.lcorr：当前步骤dzb(k)内的DGD分数

        % 对两个主态分量应用不同的相位延迟
        uux = exp(-1j * (gvd_beta + pmd_beta)) .* uux;
        uuy = exp(-1j * (gvd_beta - pmd_beta)) .* uuy;

        % 3> 通过乘以旋转矩阵R，将信号变换回原始基
        sig_fft_x = mat_rot(1, 1) * uux + mat_rot(1, 2) * uuy;
        sig_fft_y = mat_rot(2, 1) * uux + mat_rot(2, 2) * uuy;

        prop_dz = prop_dz + pmd_dz;  % 累计已处理长度
    end


else  % 仅GVD（无PMD）
    % 直接应用GVD相位旋转
    sig_fft_x = exp(-1j * phase_factor_freq * dz) .* sig_fft_x;
    sig_fft_y = exp(-1j * phase_factor_freq * dz) .* sig_fft_y;
end

% 将频谱移回标准格式
sig_fft_x = ifftshift(sig_fft_x);
sig_fft_y = ifftshift(sig_fft_y);

% 转换回时域
sig_x = ifft(sig_fft_x);
sig_y = ifft(sig_fft_y);

% 返回处理后的信号
sigout(:,1) = sig_x;
sigout(:,2) = sig_y;
end
