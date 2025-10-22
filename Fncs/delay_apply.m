function out = delay_apply(data, fs, delay)
% 对两路进行延时的判定叠加
% Apply a <delay> to an input vector <data> with samplerate <fs>.
% <delay> is given in seconds (does not have to be multiple of the
% sample interval).
% If data is complex, 对所有信号进行延迟。

% determine number of samples
n = length(data);
% make sure the input data has the right format (row-vector)
data = reshape(data, 1, n);
% the algorithm works only for even number of samples
if (mod(n,2) ~= 0)
    n = 2 * n;
    data = repmat(data, 1, 2);
    dflag = 1;
else
    dflag = 0;
end
% convert to frequency domain
fdata = fftshift(fft(real(data)))/n;
% create linear phase vector (= delay)
phd = [-n/2:n/2-1]/n*2*pi*(delay*fs);
% convert it into frequency domain
fdelay = exp(1j*(-phd));
% apply delay (convolution ~ multiplication)
fresult = fdata .* fdelay;
% ...and convert back into time domain
result = real(ifft(fftshift(fresult)))*n;


% convert to frequency domain
fdata_imag = fftshift(fft(imag(data)))/n;
% apply delay (convolution ~ multiplication)
fresult_imag = fdata_imag .* fdelay;
% ...and convert back into time domain
result_imag = real(ifft(fftshift(fresult_imag)))*n;
out = complex(result, result_imag);

if (dflag)
    out = out(1:n/2);
end
