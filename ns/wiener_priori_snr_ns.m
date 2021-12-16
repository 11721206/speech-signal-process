function wiener_priori_snr_ns(src_fn, target_fn)

if nargin<2
   fprintf('Usage: spectralsub noisyfile.wav outFile.wav \n\n');
   return;
end

[x,sr]=audioread(src_fn);

%% 首先设置参数

len = floor(20 *  sr / 1000); % 帧长 20ms
if rem(len,2)==1, len=len+1; end;
PERC=50; % 每帧的 overlap 比例, 即 10ms 帧移
lenshift = floor(PERC * len / 100);
len2 = len - lenshift;

mu = 0.98;% 噪声谱平滑因子
vad_thresh=0.15; %静音检测阈值
a_dd = 0.98; % Direct dicision 因子(直接判决因子)

%%%% 预处理
win=hamming(len); %对语音信号加汉明窗
U = (win' * win) / len; % 归一化因子

% 假设前面几帧为噪声
noise_frames = 5;
j = 1;
noise_ps = zeros(len, 1);
for i = 1:noise_frames
  noise_ps = noise_ps + abs(fft(win .* x(j: j + len - 1), len).^2)/(len * U);
  j = j + len2;
end
noise_ps = noise_ps / noise_frames;

%% 预先设定变量
k = 1;
img = sqrt(-1);
x_overlap = zeros(lenshift, 1);
nframes = floor(length(x) / len2) - 1;
xfinal = zeros(nframes * len2, 1);

%% ==================== 开始处理  ====================
% 逐帧处理
n_start = 1;
for i = 1:nframes
  frame_signal = win .* x(n_start:n_start + len - 1); % 取一帧并加窗
  spectral = fft(frame_signal, len); % 变换到频域  
  noisy_ps = abs(spectral).^2/(len * U);
  
  % 根据后验信噪比(带噪语音和噪声之间的比值)估计先验信噪比(Direct Dicision)
  if (i == 1) % 初始化后验信噪比
    posteri = noisy_ps ./ noise_ps;
	posteri_prime = posteri - 1;
	posteri_prime(find(posteri_prime < 0)) = 0;
	priori = a_dd + (1 - a_dd) * posteri_prime;
  else
    posteri = noisy_ps ./ noise_ps;
	posteri_prime = posteri - 1;
	posteri_prime(find(posteri_prime < 0)) = 0;
	priori = a_dd * (G_prev.^2) .* posteri_prev + (1 - a_dd) * posteri_prime;
  end
  
  log_sigma_k = posteri.*priori./(1 + posteri) - log(1 + priori);
  vad_dicision = sum(log_sigma_k) / len;
  if (vad_dicision < vad_thresh)
    %这个时候只有噪声,没有语音，更新噪声谱
	noise_ps = mu * noise_ps  + (1 - mu) * noisy_ps;
	vad(n_start : n_start + len -1) = 0;
  else
    vad(n_start : n_start + len -1) = 1;
  end
  
  % 更新平方根维纳滤波器系数 H = sqrt(priori / (priori + 1)) 
  G = sqrt(priori ./ (priori + 1));
  
  enhanced = ifft(spectral .* G, len);
  if (i == 1)
    xfinal(n_start : n_start + len2 - 1) = enhanced(1:len2);	
  else
    xfinal(n_start : n_start + len2 - 1) = overlap + enhanced(1:len2);
  end
  
  overlap = enhanced(len2 + 1: len);
  n_start = n_start + len2;
 
  G_prev = G;
  posteri_prev = posteri;    
end 

xfinal(n_start : n_start + len2 - 1) = overlap;

audiowrite(target_fn, xfinal, sr, 'BitsPerSample', 16); 

end