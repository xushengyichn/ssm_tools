% 定义截止频率
wc = 2*pi*10; % 10 Hz cutoff frequency

% 定义复频率范围，包括实部和虚部
[s_real, s_imag] = meshgrid(-50:0.1:50, -50:0.1:50);
s = s_real + 1j*s_imag;

% 定义传递函数
H_s = 1 ./ (1 + s/wc);

% 定义频率范围
w = -50:0.1:50;

% 定义频响函数
H_w = 1 ./ (1 + 1j*w/wc);

% 绘制传递函数和频响函数
figure;
subplot(2,1,1);
mesh(s_real, s_imag, 20*log10(abs(H_s)));
xlabel('Real Part');
ylabel('Imaginary Part');
zlabel('Magnitude (dB)');
title('Transfer Function H(s)');

subplot(2,1,2);
plot(w, 20*log10(abs(H_w)));
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Frequency Response H(jw)');
