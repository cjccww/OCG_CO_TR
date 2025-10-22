function errors=godardTED(blk,NFFT,method)
% 频域TED实现
switch (method)

    case 'Modity'
        % --- 步骤4.3: 计算Godard TED误差 ---
        error_val = sum(imag( blk(1:floor(NFFT/2)) .* conj( blk(floor(NFFT/2)+1:end) ) ) );
    case 'Conventional'
        % 两种误差计算方法
        error_val= (1/(2*pi)) * angle( sum( blk(1:NFFT/2) .* conj( blk(NFFT/2+1:end) ) ) );
end
% 存储误差
errors = error_val/NFFT;

end % function