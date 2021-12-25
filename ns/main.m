%-------------------------------------------------------------------------
%                            
%                            各种降噪示例@Matlab 2016b  
%               
%-------------------------------------------------------------------------

clc;clear all;close all;

%% 谱减法
source_fn = '../data/sp02_train_sn5.wav';
spectralsub_fn = '../data/spectralsub.wav';
spectralsub(source_fn, spectralsub_fn);

%% 基于先验信噪比估计的维纳算法
source_fn = '../data/sp02_train_sn5.wav';
wiener_priori_snr_ns_fn = '../data/wiener_priori_snr_ns.wav';
wiener_priori_snr_ns(source_fn, wiener_priori_snr_ns_fn);

%% 各种噪声估计降噪算法
% 目前支持的算法:mcra
source_fn = '../data/sp02_train_sn5.wav';
target_fn = '../data/noise_estimate_denoise.wav';
noise_estimate_denoise(source_fn, target_fn, 'mcra');
