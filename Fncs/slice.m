function [z] = slice(y, M)
% 对输入信号进行星座图切片（解调），将接收信号映射到最接近的星座点[判决函数]
% 输入:
%   y - 接收的信号（可以是实数或复数）
%   M - 星座图的阶数（符号数量）
% 输出:
%   z - 切片后的信号（最接近的星座点）

if (isreal(y))
    % 处理实数信号（如PAM调制）
    
    % 将输入信号的实部进行偏移和缩放，四舍五入得到理想星座点索引
    % 公式解释：将范围[- (M-1), M-1]映射到[0, M-1]并取整
    z_index = round( ((real(y) + (M-1)) ./ 2) );
    
    % 裁剪超出有效范围的索引值，确保索引在[0, M-1]区间内
    z_index(z_index <= -1) = 0;           % 小于0的索引设为0
    z_index(z_index > (M-1)) = M-1;       % 大于M-1的索引设为M-1
    
    % 重新生成切片后的符号，将索引映射回原始星座点范围
    z = z_index*2 - (M-1);
    
else
    % 处理复数信号（如QAM调制）
    
    % 复数星座阶数的平方根，实部和虚部分别为M_bar阶
    M_bar = sqrt(M);
    
    % 处理实部：偏移、缩放并四舍五入得到实部星座点索引
    z_index_re = round( ((real(y) + (M_bar - 1)) ./ 2) );
    
    % 处理虚部：偏移、缩放并四舍五入得到虚部星座点索引
    z_index_im = round( ((imag(y) + (M_bar - 1)) ./ 2) );

    % 裁剪实部索引，确保在有效范围[0, M_bar-1]内
    z_index_re(z_index_re <= -1)       = 0;
    z_index_re(z_index_re > (M_bar-1)) = M_bar-1;
    
    % 裁剪虚部索引，确保在有效范围[0, M_bar-1]内
    z_index_im(z_index_im <= -1)       = 0;
    z_index_im(z_index_im > (M_bar-1)) = M_bar-1;

    % 重新生成切片后的复数符号，将实部和虚部索引映射回原始星座范围
    z = (z_index_re*2 - (M_bar-1)) + 1j*(z_index_im*2 - (M_bar-1));
end
end