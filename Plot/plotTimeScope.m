function plotTimeScope(waveformBuffer,currentIndex,bufferLength,...
    mu,k,hLine,hAx)
clc;clear;close all;
% 动态显示数据示例(不适合写成函数调用的形式，更新速度太慢)

% 假设mu是一个随时间变化的信号序列，这里先生成示例数据
numPoints = 50000;  % 总数据点数
mu = 0.5 + 0.5*sin(0.01*(1:numPoints)) + 0.1*randn(1, numPoints);  % 示例信号

% 初始化图形
figure;
hFig = gcf;
hFig.Name = 'Fractional Interval';  % 窗口标题
hAx = axes(hFig);
hold on;
% 设置坐标轴范围（对应原TimeScope的配置）
xlim(hAx, [1 1e4]);  % 时间跨度1e4个点
ylim(hAx, [-2 2]);   % Y轴范围[-1,1]


legendArrary=[];
flag=struct();
flag.LegendON_OFF=0;
Plotter('Fractional Interval','Sample Index','Amplitude',[1 1e4],[-2 2],legendArrary,flag);



% 初始化波形对象 - 使用初始点避免空对象
hLine = plot(hAx, 1, mu(1), 'b');  % 先绘制第一个点

% 用于存储最近的波形数据（对应BufferLength）
bufferLength = 1e4;
waveformBuffer = zeros(1, bufferLength);
currentIndex = 0;

% 动态更新波形
for k = 1:numPoints
    % 更新缓冲区（保留最新的bufferLength个点）
    currentIndex = currentIndex + 1;
    if currentIndex > bufferLength
        % 当缓冲区满时，移除最旧的数据，添加新数据
        waveformBuffer = [waveformBuffer(2:end), mu(k)];
        xlim(hAx, [currentIndex-1e4+1 currentIndex]);
    else
        % 缓冲区未满时，直接添加新数据
        waveformBuffer(currentIndex) = mu(k);
    end
    if currentIndex < 1e4
        % 确定当前显示的范围（最近1e4个点）
        displayStart = max(1, currentIndex - 1e4 + 1);
        displayData = waveformBuffer(displayStart:currentIndex);
        displayX = displayStart:currentIndex;
    else
        displayData = waveformBuffer;
        displayX=currentIndex -1e4 + 1:currentIndex;
    end

    % 使用下标赋值方式更新数据
    hLine.XData(1:length(displayX)) = displayX;
    hLine.YData(1:length(displayX)) = displayData;

    % 刷新图形
    drawnow limitrate;  % 限制刷新速率，提高性能
end
