function spectralsub(src_fn, target_fn)

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

thres = 3;
alpha=2.0; 
FLOOR=0.002;
G=0.9;

%%%% 预处理
win=hamming(len); %对语音信号加汉明窗
winGain = len2 / sum(win); % 

nFFT=2 * 2^nextpow2(len);
noise_mean = zeros(nFFT, 1); % 
% 假设最开始的 5 帧是静音段
j = 1;
for i = 1:5
  noise_mean = noise_mean + abs(fft(win .* x(j: j + len - 1), nFFT));
  j = j + len;
end
noise_mean = noise_mean / 5;

%% 预先设定变量
k = 1;
img = sqrt(-1);
x_overlap = zeros(lenshift, 1);
nframes = floor(length(x) / len2) - 1;
xfinal = zeros(nframes * len2, 1);

%% ==================== 开始处理  ====================
% 逐帧处理
for i = 1:nframes
  frame_signal = win .* x(k:k + len - 1); % 取一帧并加窗
  spectral = fft(frame_signal, nFFT); % 变换到频域
  mag = abs(spectral);
  theta = angle(spectral);
  
  snr = 10 * log10(norm(mag,2)^2 / norm(noise_mean, 2)^2);%计算信噪比
  
  if alpha==1.0  % 幅度
      beta=berouti1(snr);
  else % 功率
      beta=berouti(snr);
  end
  
  %根据公式计算 D(w) = 带噪信号谱(功率或者幅度)- beta * 噪声谱(幅度或者功率)
  d = spectral.^alpha - beta * noise_mean.^alpha; % 得到纯净语音
  diffw = d - FLOOR * noise_mean.^alpha; % 相减后的如果大于某个门限，那就是不做处理，如果小于该门限，之前是暴力设置为0，此处将噪声谱作为该门限值,减少音乐噪声带来的影响
  z=find(diffw <0);  
  if ~isempty(z)
    d(z)=FLOOR*noise_mean(z).^alpha;
  end
  
  % ---------- 更新噪声谱 ----------
  if snr < thres
    noise_temp = G * noise_mean.^alpha + (1 - G) * mag.^alpha;
	noise_mean = noise_temp.^(1/alpha);
  end
  
  d(nFFT/2+2 : nFFT) = flipud(d(2:nFFT/2)); % 取另外一半的频谱
  x_d = d.^(1/alpha) .*(cos(theta) + img * sin(theta)); % 得到 x = 幅度谱 * 相位谱
  xi = real(ifft(x_d)); % 逆变换
  xfinal(k : k + len2 - 1) = x_overlap + xi(1:lenshift); % overlap
  x_overlap = xi(1+lenshift : len); % 把overlap的语音放在x_overlap中，供下一帧语音进行计算
  
  k = k + len2;
end 

audiowrite(target_fn, xfinal * winGain, sr, 'BitsPerSample', 16); % 去窗并写文件


%--------------------------------------------------------------------------

function a=berouti1(SNR)

if SNR>=-5.0 & SNR<=20
   a=3-SNR*2/20;
else
   
  if SNR<-5.0
   a=4;
  end

  if SNR>20
    a=1;
  end
  
end

function a=berouti(SNR)

if SNR>=-5.0 & SNR<=20
   a=4-SNR*3/20; 
else
   
  if SNR<-5.0
   a=5;
  end

  if SNR>20
    a=1;
  end
  
end
