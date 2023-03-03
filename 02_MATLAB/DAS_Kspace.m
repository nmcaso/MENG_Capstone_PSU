            clc; close all;

addpath("03_Classes\");
x_arr       = -0.1:0.001:0.01;
delay_arr   = 0.4*(1:length(x_arr));

freq        = 100;
tau         = 2*pi;
fs          = 1000;
N           = 100;
time        = (0:N-1)./fs;

s_tx        = cos(tau*freq*time'+delay_arr);
s_tk        = fft(s_tx,[],2);
s_omgx      = fft(s_tx);
s_omgk      = fftn(s_tx);

figure; tiledlayout(1,4); nexttile;
imagesc(20*log10(abs(s_tk))); colormap('pucolors.cvidis'); daspect([1 1 1]); nexttile;
imagesc(20*log10(abs(s_omgx))); colormap('pucolors.cvidis'); daspect([1 1 1]); nexttile;
imagesc(20*log10(abs(s_omgk))); colormap('pucolors.cvidis'); daspect([1 1 1]); nexttile;
imagesc(s_tx); colormap(pucolors.cvidis); daspect([1 1 1]);