function [slope_s,zero_x] = zero_slope(phase, s_curve)
    % 方法1：找到过零点附近的区域进行线性拟合
    zero_threshold = 0.01;  % 定义"接近零点"的阈值
    near_zero_indices = find(abs(s_curve) < zero_threshold);
    
    if length(near_zero_indices) < 10
        % 方法2：如果直接找零点附近点不够，使用相位接近0的区域
        phase_threshold = 0.05;  % 使用相位接近0的区域
        near_zero_indices = find(abs(phase) < phase_threshold);
    end
    
    % 对零点附近区域进行线性拟合
    p = polyfit(phase(near_zero_indices), s_curve(near_zero_indices), 1);
    slope_s = p(1);  % 斜率就是S曲线在零点的增益
    zero_x = -p(2)/p(1);
    
    % 可视化验证
%     figure;
%     plot(phase, s_curve, 'b-', 'LineWidth', 2); hold on;
%     plot(phase(near_zero_indices), s_curve(near_zero_indices), 'ro', 'MarkerSize', 4);
%     plot(phase, polyval(p, phase), 'r--', 'LineWidth', 1.5);
%     xlabel('相位误差 (UI)'); ylabel('鉴相器输出');
%     title('S曲线及零点斜率提取');
%     legend('S曲线', '用于拟合的点', '线性拟合', 'Location', 'best');
%     grid on;
end