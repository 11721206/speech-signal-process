%-------------------------------------------------------------------------
%                            
%                            各种降噪示例@Matlab 2016b  
%               
%-------------------------------------------------------------------------

clc;clear all;close all;

%% 谱减法
source_fn = 'sp02_train_sn5.wav';
spectralsub_fn = 'spectralsub.wav';
spectralsub(source_fn, spectralsub_fn);
