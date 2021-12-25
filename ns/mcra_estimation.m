function parameters = mcra_estimation(ns_ps, parameters)

as = parameters.as;
ad = parameters.ad;
ap = parameters.ap;
pk = parameters.pk;
delta = parameters.delta;
L = parameters.L;
n = parameters.n;
len = parameters.len;
noise_ps = parameters.noise_ps;
P = parameters.P;
Pmin = parameters.Pmin;
Ptmp = parameters.Ptmp;

% 计算局部能量
P = as * P + (1 - as) * ns_ps;

% 计算最小能量
if rem(n, L) == 0
  Pmin = min(Ptmp, P);
  Ptmp = P;
else
  Pmin = min(Pmin, P);
  Ptmp = min(Ptmp, P);
end

%计算当前能量和做小局部能量的比值
Srk = P./Pmin;

% 平滑语音存在概率 
Ikl = zeros(len,1);
Ikl(find(Srk > delta)) = 1; 
pk = ap * pk + (1 - ap) * Ikl;

% 计算论文中的公式 5
adk = ad + (1 - ad) * pk;

% 更新噪声谱 公式 4
noise_ps = adk.*noise_ps + (1 - adk).*ns_ps;

% 记录当前帧的一些信息
parameters.pk = pk;
parameters.n = n+1;
parameters.noise_ps = noise_ps;
parameters.P = P;
parameters.Pmin = Pmin;
parameters.Ptmp = Ptmp;

end