function [y,phEst] = symbolTimingPLL_Feedforward(signal,sps,szBlock,bias,ted_type,intpl,const,debug_s, debug_r,debug_r_fast)
% Feedforward timing phase estimation and recovery. The input
% parameters are defined as follows:
%
%   szBlock: block size
%   estMeth: estimator method
%   intMeth: interpolation method
%   decFlag: decimate flag

% 参数校验与默认值设置
if (nargin < 10)
    debug_r_fast = 0; % 默认关闭实时调试示波器
end
if (nargin < 9)
    debug_r = 0; % 默认关闭实时调试示波器
end
if (nargin < 8)
    debug_s = 0; % 默认关闭静态调试绘图
end
if nargin<6
    intpl = 1;
end
if nargin<5
    ted_type = 'lee';
end
if nargin<4
    bias = 1.0;
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

%%
% 输入信号
x = signal;


% get the size of input
nSamples = size(x,1);
kk = size(x,2);

% 余数
temp = mod(nSamples,szBlock);
if temp
    x = [ x; zeros(szBlock-temp,kk)];
end
% Feedforward PLL
[y,phEst,mu,e]= PLL(x,sps,bias,szBlock,intpl,ted_type,const,debug_r_fast);


%% 算法运行图例及参数变化图
if (debug_s)
    figure
    plot(e(:,1))
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
    title('TED output')

    figure
    plot(phEst)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(mu(:,1), '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
end


%%

function [symRtm,tau,nu,e]= PLL(x,sps,bias,szBlk,intpl,estMeth,const,debug_r_fast)
% Foreford_PLL This is NOT a phase lock loop
% 快速调试图
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
% Initialize 

% 初始防止延迟跳过超过一个/多个符号
temp = 0;
while (temp+1/sps-1/2)*(temp-1/2) >= 0
    px = x(1:szBlk,1);
    if sum(strcmpi(estMeth,{'sln','fln','avn'})) && sps == 2
        % interpolation to get 4 samples per symbol using interpft
        px = interpft(px, 2*length(px));
    end
    temp = TED(px, estMeth, bias);
    temp = -1/2/pi * angle(temp);
    id = floor((0.5-temp)*sps);
    x(:,1) = circshift(x(:,1), -id);
end

temp = 0;
while (temp+1/sps-1/2)*(temp-1/2) >= 0
    px = x(1:szBlk,2);
    if sum(strcmpi(estMeth,{'sln','fln','avn'})) && sps == 2
        % interpolation to get 4 samples per symbol using interpft
        px = interpft(px, 2*length(px));
    end
    temp = TED(px, estMeth, bias);
    temp = -1/2/pi * angle(temp);
    id = floor((0.5-temp)*sps);
    x(:,2) = circshift(x(:,2), -id);
end

% x = [x;zeros(sps*10,2)];
pointer = [1 1]; % first index of basepoint set
k = 1;
ind = 0:szBlk-1;
while max(pointer)+szBlk+sps-2 <= size(x,1)
    
    px1 = x(pointer(1)+ind,1);
    px2 = x(pointer(2)+ind,2);
    if sum(strcmpi(estMeth,{'sln','fln','avn'})) && sps == 2
        % interpolation to get 4 samples per symbol using interpft
        px1 = interpft(px1, 2*length(px1));
        px2 = interpft(px2, 2*length(px2));
    end
    
    % ================= TED =================== %
    sumcomp(k,1) = TED(px1, estMeth, bias);
    sumcomp(k,2) = TED(px2, estMeth, bias);
    
    % ================= Loop Filter =================== %
    % % smooth the complex signal
    summation(k,1) = loopfilter(sumcomp(:,1));
    summation(k,2) = loopfilter(sumcomp(:,2));
    tau(k,1) = -1/2/pi*angle(summation(k,1));
    tau(k,2) = -1/2/pi*angle(summation(k,2));
    
    % TED error 
    e(k,1)=tau(k,1);
    e(k,2)=tau(k,2);

    % ================= Control =================== %
    if k == 1 
        tau1 = 0;
    else 
        tau1 = tau(k-1,1);
    end
    tau2 = tau(k,1);
    [pointer(1), nu(k,1), tau(k,1)] = control(tau1, tau2, pointer(1), sps);
    
    if k == 1 
        tau1 = 0;
    else 
        tau1 = tau(k-1,2);
    end
    tau2 = tau(k,2);
    [pointer(2), nu(k,2), tau(k,2)] = control(tau1, tau2, pointer(2), sps);
    
    % ================= Interpolation =================== %
    % Sampling instants for tau calculation
    for n = 1:szBlk
        interpolantee = x(pointer(1)+n-1:pointer(1)+n+sps-2,1); % 可以 每次 通用地取出 sps 个连续采样点
        datInt1(n,k) = interpolate(interpolantee,nu(k,1),intpl);
        interpolantee = x(pointer(2)+n-1:pointer(2)+n+sps-2,2);
        datInt2(n,k) = interpolate(interpolantee,nu(k,2),intpl);
    end

    if (debug_r_fast)
        step(hScope, datInt1(:,k)) % 更新星座图
        step(hTScopeCounter, nu(k,1)); % 更新μ值曲线
    end
    % 基点
    pointer = pointer+szBlk;
    k = k+1;
end

% format output
symRtm = [datInt1(:), datInt2(:)];


%% 插值函数
function y = interpolate(x,nu,intpl)
% Interpolate signal samples
sps = size(x,1);
switch (intpl)
    case 0
        y = x(1);
    case 1
        if sps == 2
            y = x(1,:) + nu*(x(2,:)-x(1,:));
        else
            error('Lee only need 2sps');
        end
    case 2
        if sps ==4 
        a = 0.5;
        b = [0 1 0 0;...
            -a a-1 a+1 -a;...
            a -a -a a]; % b(i,l)
        v = b*x;
        y = (v(3)*nu+v(2))*nu+v(1);
        else
            error('SLN/FLN only need 4sps');
        end
    case 3
        if sps ==4 
        b = [0 1 0 0;...
            -1/3 -1/2 1 -1/6;...
            1/2 -1 1/2 0;...
            -1/6 1/2 -1/2 1/6]; % b(i,l)
        v = b*x;
        y = ((v(4)*nu+v(3))*nu+v(2))*nu+v(1);
        else
            error('SLN/FLN only need 4sps');
        end
end


%% 环路滤波
function y = loopfilter(x)
y = x(end);


%% NCO 控制
function [pointer, nu, tau2] = control(tau1, tau2, pointer, sps)
if floor((0.5-tau2)*sps)
    if tau1-tau2>0.75
        pointer = pointer+1;    % skip one sample
        tau2 = tau2+(sps-1)/sps;
        nu = tau2*sps;
    elseif tau1-tau2<=0.25
        pointer = pointer-1;     % wait for one sample
        tau2 = tau2+1/sps;
        nu = tau2*sps;
    elseif tau1-tau2<=0.75 && tau1-tau2>0.5
        pointer = pointer+2;     % skip 2 sample
        tau2 = tau2+(sps-2)/sps;
        nu = tau2*sps;
    elseif tau1-tau2<=0.5 && tau1-tau2>0.25
        pointer = pointer-2;     % wait for 2 sample
        tau2 = tau2+2/sps;
        nu = tau2*sps;
    end
else
    nu = tau2*sps;
end
nu = mod(nu,1);



%% TED Output
function s = TED(px, TED, g)

if nargin<3
    g = 1;
end

switch TED
    case 'fln'
        s = ted_fln(px);
    case 'sln'
        s = ted_sln(px);
    case 'avn'
        s = ted_avn(px);
    case 'lee'
        s = ted_lee(px,g);
    case 'godard'
        s = ted_godard(px);        
end

%% ted type

function s = ted_fln(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px).^4 .* ex.' );

function s = ted_sln(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px).^2 .* ex.' );

function s = ted_avn(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px) .* ex.' );

function s = ted_lee(px,g)
L = length(px);
n = 1:L;
% cosine part
ex1 = (-1).^(n-1);          % ex1 = exp(-1j.*(ii-1).*pi);
sum_1 = sum( abs(px).^2 .* ex1.' );
% sine part
ex2 = 1j * (-1).^(n-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
xh = px(2:end);
xx = px(1:end-1);
ex2 = ex2(1:end-1);
sum_2 = sum( real(conj(xx).*xh) .* ex2.' );
% with biasing
s = g*sum_1 + sum_2;
s = conj(s);

function s = ted_godard(px)
N = length(px);
X = fft(px);
Z = xcorr(X);
s = Z(N/2);
