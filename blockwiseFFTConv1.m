function y = blockwiseFFTConv1(x, h, NFFT, freqDomainFilter)
% BLOCKWISEFFTCONV Blockwise convolution in the frequency domain using the overlap-and-save FFT method.
%
% Parameters:
%   x : Input signal (vector)
%   h : Filter impulse response or frequency response (vector)
%   NFFT : FFT size to be used. Must be greater than the length of the filter.
%          If empty, it will be set to the next power of 2 greater than or equal
%          to the length of the filter. Default is [].
%   freqDomainFilter : If true, `h` is assumed to be the frequency response of the filter.
%                      If false, the FFT of `h` will be computed. Default is false.
%
% Returns:
%   y : The filtered output signal.

    if nargin < 4
        freqDomainFilter = false;
    end
    if nargin < 3
        NFFT = [];
    end

    sigLen = length(x);  % length of the input signal
    M = length(h);      % length of the filter impulse response
    D = floor((M - 1) / 2);  % filter delay

    % Determine NFFT if not provided
    if isempty(NFFT)
        NFFT = 2^ceil(log2(M));
    end

    % Validate NFFT
    if NFFT < M
        error('FFT size is smaller than filter length');
    end

    L = NFFT - M + 1;  % block length required

    % Preprocess filter
    if freqDomainFilter
        % Convert frequency response to time domain, apply fftshift, then pad
        h_time = ifft(h);
        h_time_shifted = fftshift(h_time);
        h_padded = [h_time_shifted; zeros(L - 1, 1)];
    else
        % Pad filter with zeros
        h_padded = [h(:); zeros(L - 1, 1)];
    end

    H = fft(h_padded);  % frequency response

    discard = M - 1;  % number of samples to be discarded after IFFT
    numBlocks = ceil(sigLen / L);  % total number of FFT blocks
    padLen = numBlocks * L - sigLen;  % pad length for complete blocks

    % Pad signal with zeros (end padding + delay compensation)
    x_padded = [x(:); zeros(padLen + D, 1)];

    % Pre-allocate output
    y = zeros(length(x_padded), 1);

    % Add beginning padding for overlap-and-save
    x_padded = [zeros(M - 1, 1); x_padded];

    % Overlap-and-save blockwise processing
    start_idx = 1;
    end_idx = NFFT;

    for blk = 1:numBlocks
        % Extract current block and compute FFT
        x_block = x_padded(start_idx:end_idx);
        X = fft(x_block);

        % Multiply in frequency domain and compute IFFT
        y_blk = ifft(X .* H);

        % Store valid part (discard overlap samples)
        output_start = (blk - 1) * L + 1;
        output_end = blk * L;
        y(output_start:output_end) = y_blk(discard + 1:end);

        % Move to next block
        start_idx = start_idx + L;
        end_idx = start_idx + NFFT - 1;
    end

    % Remove padding and return final result
    y_final = y(D + 1:end - padLen);

    % Return real part if input was real
    if isreal(x)
        y_final = real(y_final);
    end
end