%-------------------------------------------------------------------------
%                            
%                            波束响应图@Matlab 2016b  
%                            encoding : UTF-8
%-------------------------------------------------------------------------

% 均匀线阵; 阵元个数 N = 11, 阵元间距 半波长, 声速 c = 340, 均匀加权1/N 频率 1700Hz

%-------------------------------------------------------------------------
%                            第一部分
%-------------------------------------------------------------------------

%% 画出频率波束响应图; 主波束方向为 90
close all; clc; clear all;
N = 11;
c = 340;
f = 1700;
lamda = c/f;  
d = lamda / 2;
xPos = -0.25:0.05:25;
yPos = zeros(1, numel(xPos));
zPos = zeros(1, numel(xPos));
elementWeights = ones(1, numel(xPos)) / N;

% 线性扫描的角度范围
thetascanAngle = 0:0.001:3.1415;  % 从 0 到 pi
phiScanAngles = 0;

B = zeros(1, numel(thetascanAngle));

for i = 1 : numel(thetascanAngle)
    p = exp(-1j * (N - 1) / 2 * 2 * pi * d / lamda * cos(thetascanAngle(i)));
	for k = 1 : N
	    B(i) = B(i) + elementWeights(k) * exp(1j * (k - 1) * 2 * pi * d / lamda * cos(thetascanAngle(i)));
	end
	B(i) = B(i) * p;
end
hpol = polardb(thetascanAngle, 10 * log10(abs(B).^2), -50);title('ULA polar beam pattern');
figure;plot(thetascanAngle, 10 * log10(abs(B).^2));title('ULA rect beam pattern'); grid on;

%-------------------------------------------------------------------------
%                            第二部分
%-------------------------------------------------------------------------
%% 主波束调向至 120 度

B1 = zeros(1, numel(thetascanAngle));
target_theta = 2*pi/3;
for i = 1 : numel(thetascanAngle)
    p = exp(-1j * (N - 1) / 2 * 2 * pi * d / lamda * (cos(thetascanAngle(i)) - cos(target_theta)));
	for k = 1 : N
	    B1(i) = B1(i) + elementWeights(k) * exp(1j * (k - 1) * 2 * pi * d / lamda * (cos(thetascanAngle(i)) - cos(target_theta)));
	end
	B1(i) = B1(i) * p;
end
figure;hpol = polardb(thetascanAngle, 10 * log10(abs(B1).^2), -50);title('ULA polar beam pattern-主波束调向');
figure;plot(thetascanAngle, 10 * log10(abs(B1).^2));title('ULA rect beam pattern-主波束调向'); grid on;

%-------------------------------------------------------------------------
%                            第三部分
%-------------------------------------------------------------------------
%% 零点调向 在 135 度形成零点调向
null_angle = 3*pi/4;
