function [Eo, t_nco_values] = gardnerClockRecovery(Ei, ki,kp,intpl,isNyquist,lpad,maxPPM,debug_s)
% GARDNERCLOCKRECOVERY Perform clock recovery using Gardner's algorithm with a loop PI filter
%
% Inputs:
%   Ei     - Input array representing the received signal
%   param  - Structure containing clock recovery parameters (optional)
%     .kp          - Proportional gain for the loop filter (default: 1e-3)
%     .ki          - Integral gain for the loop filter (default: 1e-6)
%     .isNyquist   - Is the pulse shape a Nyquist pulse? (default: true)
%     .returnTiming- Return estimated timing values (default: false)
%     .lpad        - Length of zero padding at the end of input vector (default: 1)
%     .maxPPM      - Maximum clock rate expected deviation in PPM (default: 500)
%
% Outputs:
%   Eo     - Recovered signal
%   t_nco_values - Timing values (optional, only if returnTiming is true)

% Set default parameters if not provided
if (nargin < 4)
    intpl = 3;
end
if (nargin < 5)
    isNyquist = 0;
end
if (nargin < 6)
    lpad = 1;
end
if (nargin < 7)
    maxPPM = 500;
end
if (nargin < 8)
    debug_s = 0;
end

% Ensure Ei is a column vector if 1D
if isvector(Ei)
    Ei = Ei(:);
end

% Pad the input signal
[nSamples, nModes] = size(Ei);
Ei = [Ei; zeros(lpad, nModes)];

% Initialize output vector accounting for maximum clock deviation
Ln = floor((1 - maxPPM / 1e6) * (nSamples+lpad));
Eo = zeros(Ln, nModes);
t_nco_values = zeros(Ln, nModes);

last_n = 0;
fprintf('Running clock recovery...\n');
% 赋值前两个元素(该写法并不需要对前两个数值进行赋值)
%Eo(1:2,:)=Ei(1:2,:);
for indMode = 1:nModes
    intPart = 0;
    t_nco = 0;
    %mu = 0;
    % n为输出数组的索引，m为输入数组的索引
    n = 3;
    m = 3;

    % loopFilterOut=1;
    while n < Ln  && m < size(Ei, 1) - 1
        % Interpolator function call
        Eo(n, indMode) = interpolator(Ei(m-2:m+1, indMode), t_nco,intpl);

        if mod(n, 2) == 0  % Even sample index
            if isNyquist
                ted(n,indMode) = gardnerTEDnyquist(Eo(n-2:n, indMode));
            else
                ted(n,indMode) = gardnerTED(Eo(n-2:n, indMode));
            end

            % Loop PI Filter
            intPart = ki * ted(n,indMode) + intPart;
            propPart = kp * ted(n,indMode);
            loopFilterOut(n,indMode) = propPart + intPart;

            t_nco = t_nco - loopFilterOut(n,indMode);
        end

        % NCO clock gap adjustment
        if t_nco > 1
            t_nco = t_nco - 1;
            n = n - 1;
        elseif t_nco < -1
            t_nco = t_nco + 1;
            n = n + 2;
            m = m + 1;
        else
            n = n + 1;
            m = m + 1;
        end

        t_nco_values(n, indMode) = t_nco;

    end

    if n > last_n
        last_n = n;
    end

end

% Trim output to actual processed length
Eo = Eo(1:last_n, :);
t_nco_values = t_nco_values(1:last_n, :);


% 算法运行图例及参数变化图
if (debug_s)
    % Calue ppm
    [ppm,peak_locs] = calcClockDrift(t_nco_values);
    % Iteration of SCO
    t_nco_convergence=[];
    for i =2:length(peak_locs)-1
        y_linear=t_nco_values(peak_locs(i)+1:peak_locs(i+1));
        p_linear = polyfit((1:length(y_linear)), y_linear, 1);  % 拟合1次多项式（线性）
        trend_fit_linear = polyval(p_linear, (1:length(y_linear)));  % 计算拟合的趋势
        y_detrend_linear_builtin = y_linear - trend_fit_linear.';  % 去趋势后的数据
        t_nco_convergence=[t_nco_convergence;y_detrend_linear_builtin];
    end
    t_nco_values1=[t_nco_values(1:peak_locs(2));t_nco_convergence+1];
    % 计算SCO
    SCO=t_nco_values1*ppm;

    error=t_nco_values(peak_locs(2)+1:peak_locs(end))-t_nco_convergence;
    RMSE = mean((error).^2);
    fprintf("RMSE  %.2f \n",RMSE)

    figure
    plot(ted(:,1))
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
    title('TED output')

    figure
    plot(loopFilterOut(:,1))
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(t_nco_values(:,1), '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')

    figure;box on;
    plot(SCO)
    title('Iteration')
    ylabel('SCO(ppm)')
    xlabel('Symbol $k$', 'Interpreter', 'latex')



end



end

%% TED Type
function ted = gardnerTED(x)
% 计算Gardner定时误差检测器的输出值
% 输入x: 大小为3的数组，表示接收信号的一段
% 输出ted: Gardner定时误差检测值

% 在MATLAB中，conj()用于计算复数共轭，real()取实部
ted = real(conj(x(2)) * (x(3) - x(1)));
end
function ted = gardnerTEDnyquist(x)
% 适用于奈奎斯特脉冲的改进型Gardner定时误差检测器
% 输入x: 大小为3的数组，表示接收信号的一段
% 输出ted: 改进型Gardner定时误差检测值

% abs()计算复数的模，.^2表示元素-wise平方运算
ted = (abs(x(2))^2) * (abs(x(1))^2 - abs(x(3))^2);
end

%% 插值函数
function y = interpolator(x,nu,intpl)
% Interpolate signal samples
sps = size(x,1);
switch (intpl)
    case 0
        y = x(1);
    case 1
        if sps == 2
            y = x(1,:) + nu*(x(2,:)-x(1,:));
        else
            error('The interpolator methord need 2sps');
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
            error('The interpolator methord need 4sps');
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
            error('The interpolator methord need 4sps');
        end
end
end