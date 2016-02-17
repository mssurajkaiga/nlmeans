data = ones(100, 100, 3);
data = imnoise(data, 'gaussian', 0, 0.3);
imwrite(data, 'input_gaussian_zeros_0_0-1.jpeg');
hdrwrite(data, 'input_gaussian_ones_0_0-3.hdr');