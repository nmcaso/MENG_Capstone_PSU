%%
% This script will run DAS and DMAS over the five data sets from MB1-MB5
%%
clear

load RCV_Blood_MB1_400_500.mat
%%
sound_speed = 1380;

wL = sound_speed/Trans.frequency/1e6

time_vec = double(PARcvData(17:end,:,1)');

sens_vec = Trans.ElementPos(:,1)*1e-3;

x_start = 70*wL;
x_end = 160*wL;
x_length = x_end - x_start;

x_vec = [x_start:1/60*1e-3:x_end]';

% y_start = Trans.ElementPos(1,1)/2/1e3;
% y_end = Trans.ElementPos(end,1)/2/1e3;
y_start = -9.6/1e3;
y_end = 9.6/1e3;
y_length = y_end - y_start;
y_vec = [y_start:1/60*1e-3:y_end];

fs = 4*5.208e6;
dt = 1/fs;

pulse_index = 0;
sens_count = length(sens_vec);
t_max = size(time_vec,2); %maximum index number for the time delay matrix
total_pixels = length(x_vec)*length(y_vec);

%%
"DAS reconstruction"
tic
DAS_image = reshape(DAS_vectorized(x_vec,y_vec,time_vec,sens_vec,sound_speed,dt,pulse_index,sens_count,t_max,total_pixels),length(x_vec),length(y_vec));
toc
%%
figure
imagesc(y_vec*1e3,x_vec*1e3,flip(abs(hilbert(DAS_image(end:-1:1,:))),1))
title('DAS hot filtered')
% axis equal
% ylim([0.022 0.036])
colormap(hot)
colorbar
pbaspect([y_length/x_length 1 1])
xlabel('Lateral [mm]'); ylabel('Depth [mm]')
caxis([0 7]*1e4)

"DAS prefiltered"
figure
imagesc(y_vec*1e3,x_vec*1e3,flip(DAS_image(end:-1:1,:),1))
title('DAS hot unfiltered')
% axis equal
% ylim([0.022 0.036])
colormap(hot)
colorbar
pbaspect([y_length/x_length 1 1])
xlabel('Lateral [mm]'); ylabel('Depth [mm]')
caxis([-6 6]*1e4)

%%
"DMAS reconstruction"
tic
DMAS_image = reshape(DMAS_vectorized(x_vec,y_vec,time_vec,sens_vec,sound_speed,dt,pulse_index,sens_count,t_max,total_pixels),length(x_vec),length(y_vec));
toc
%%
figure
imagesc(y_vec*1e3,x_vec*1e3,flip(abs(hilbert(DMAS_image(end:-1:1,:))),1))
title('DMAS hot filtered')
% ylim([0.022 0.036])
colormap(hot)
colorbar

figure
imagesc(y_vec*1e3,x_vec*1e3,flip(DMAS_image(end:-1:1,:),1))
title('DMAS hot unfiltered')
% ylim([0.022 0.036])
colormap(hot)
colorbar
%%

% "DAS movie"
% tic
% writerObj = VideoWriter('DAS_movie.avi');
% writerObj.FrameRate = 30;
% for i=1:500
%     time_vec = double(PARcvData(17:end,:,i)');
%     t_max = size(time_vec,2);
%     DAS_image = reshape(DAS_vectorized(x_vec,y_vec,time_vec,sens_vec,sound_speed,dt,pulse_index,sens_count,t_max,total_pixels),length(x_vec),length(y_vec));
%     figure(2)
%     imagesc(x_vec,y_vec,abs(hilbert(DAS_image(end:-1:1,:)')))
%     axis equal
%     colormap(jet)
%     colorbar
%     F(i) = getframe;
% end
% 
% open(writerObj)
% for i=1:length(F)
%     frame = F(i);
%     writeVideo(writerObj, frame);
% end
% close(writerObj)
% toc