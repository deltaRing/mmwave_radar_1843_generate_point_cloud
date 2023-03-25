clear all;close all;clc;
GHz = 10e9;
MHz = 10e6;
us = 10e-6;
%% 雷达参数（使用mmWave Studio默认参数）
c=3.0e8;  
B=4000*10e6;       %调频带宽
K = 68.992*10e12;  %调频斜率
T=B/K;         %采样时间
Tc=65e-6;     %chirp总周期
fs=10*10e6;       %采样率
f0=77e9;       %初始频率
lambda=c/f0;   %雷达信号波长
d=lambda/2;    %天线阵列间距
n_samples=256; %采样点数/脉冲
N=256;         %距离向FFT点数
n_chirps=64;   %每帧脉冲数
M=256;         %多普勒向FFT点数
n_RX=4;        %RX天线通道数
Q = 512;       %角度FFT
tx = 2;        %发射天线数目
rx = 4;        %接受天线数目
fnumber = 256;
PRI = 40e-3;

range_axis = linspace(0, fs * c / 2 / K, N);
velo_axis = linspace(-lambda / 4 / Tc, lambda / 4 / Tc, M);
an_axis_az = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
an_axis_el = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
fname = "2023-03-05-16-33-52.bin";

fid = fopen(fname,'rb'); 
sdata = fread(fid,'int16');    
for xx = 1:fnumber   
    %16bits，复数形式(I/Q两路)，4RX,3TX,有符号16bit，小端模式
    sdata2 = sdata((xx-1)*n_samples*n_chirps*rx*tx*2+1:xx*n_samples*n_chirps*rx*tx*2);
    %% 1843+DCA1000
    % 阵列排布
    %       rx1 rx2  rx3 rx4 
    %           =======> 方位向
    % tx1 ： []  []  []  []   ||
    % tx2 ： []  []  []  []   ||
    % tx3 ： []  []  []  []   ||
    %               俯仰向    \/
    fileSize = size(sdata2, 1);
    lvds_data = zeros(1, fileSize/2);
    count = 1;
    for i=1:4:fileSize-5
       lvds_data(1,count) = sdata2(i) + 1i*sdata2(i+2); 
       lvds_data(1,count+1) = sdata2(i+1)+1i*sdata2(i+3); 
       count = count + 2;
    end
    lvds_data = reshape(lvds_data, tx, n_samples*n_RX, n_chirps);
    
    tx_data = [];
    tx_select_rx = 2;
    rx_select_tx = 2;
    for ttxx = 1:tx
        lvds_data_tx = squeeze(lvds_data(ttxx,:,:));
        lvds_data_tx = lvds_data_tx.';
        cdata = zeros(n_RX,n_chirps*n_samples);
        for row = 1:n_RX
              for i = 1: n_chirps
                  cdata(row,(i-1)*n_samples+1:i*n_samples) = lvds_data_tx(i,(row-1)*n_samples+1:row*n_samples);
              end
              if tx_select_rx == row
                  temp_data = reshape(cdata(row,:),n_samples,n_chirps);
%                   temp_data = temp_data(:,2:end) - temp_data(:,1:end-1);
                  tx_data(:,:,ttxx) = temp_data;
              end
        end
        
        if ttxx == rx_select_tx
            data_radar_1 = reshape(cdata(1,:),n_samples,n_chirps);   %RX1
            data_radar_2 = reshape(cdata(2,:),n_samples,n_chirps);   %RX2
            data_radar_3 = reshape(cdata(3,:),n_samples,n_chirps);   %RX3
            data_radar_4 = reshape(cdata(4,:),n_samples,n_chirps);   %RX4
        end
    end
    
%     figure(6)
%     subplot(221)
%     plot(real(squeeze(data_radar_1(:,1))))
%     subplot(222)
%     plot(real(squeeze(data_radar_2(:,1))))
%     subplot(223)
%     plot(real(squeeze(data_radar_3(:,1))))
%     subplot(224)
%     plot(real(squeeze(data_radar_4(:,1))))
    
    data_radar_1 = data_radar_1(:,7:end) - data_radar_1(:,1:end-6);
    data_radar_2 = data_radar_2(:,7:end) - data_radar_2(:,1:end-6);
    data_radar_3 = data_radar_3(:,7:end) - data_radar_3(:,1:end-6);
    data_radar_4 = data_radar_4(:,7:end) - data_radar_4(:,1:end-6);
%     data_radar_1 = data_radar_1 - mean(data_radar_1, 2);
%     data_radar_2 = data_radar_2 - mean(data_radar_2, 2);
%     data_radar_3 = data_radar_3 - mean(data_radar_3, 2);
%     data_radar_4 = data_radar_4 - mean(data_radar_4, 2);
    data_radar=[];            
    data_radar(:,:,1)=data_radar_1;     %三维雷达回波数据
    data_radar(:,:,2)=data_radar_2;
    data_radar(:,:,3)=data_radar_3;
    data_radar(:,:,4)=data_radar_4;
    %% 3维FFT处理
    %距离FFT
    range_win = hamming(n_samples);   %加海明窗
    doppler_win = hamming(n_chirps);
    range_profile = [];
    for k=1:n_RX
        range_profile(:,:,k) = flip(fft(data_radar(:,:,k), N));
    end
    %多普勒FFT
    speed_profile = [];
    for k=1:n_RX
        speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    
    range_profile1 = flip(fft(data_radar_1, N));
    figure(10001)
    imagesc(abs(squeeze(range_profile(:,:,1))))
% %     speed_profile1 = fft(range_profile1, M, 2);
%     figure(10002)
%     imagesc(abs(speed_profile(:,:,1)))
    
 %     CACFAR
    detect_profile = [];
    detect_res = {}; % range, doppler, val
    [detect_profile, detect_res] = CFAR(squeeze(speed_profile(:,:,1)));
    
%     figure(10003)
%     mesh(abs(squeeze(speed_profile(:,:,1))))
%     figure(10004)
%     mesh(abs(detect_profile))


    az_rmap = squeeze(sum(abs(fft(speed_profile, Q, 3)), 2));
    angle_detect = AnglePeaksDetect(az_rmap);
    figure(10005)
    mesh(abs(az_rmap));
    
    angles = an_axis_az(angle_detect(:,2));
    ranges = range_axis(angle_detect(:,1));
    
    xx = ranges .* sin(angles);
    yy = ranges .* cos(angles);
    
    figure(10006)
    scatter(xx, yy)
    hold on
    axis([-15 15 0 20])
    
%     range_win = hamming(n_samples);   %加海明窗
%     doppler_win = hamming(n_chirps);
%     range_profile_tx = [];
%     for k=1:tx
%         range_profile_tx(:,:,k) = fftshift(fft(tx_data(:,:,k), N));
%     end
%     %多普勒FFT
%     speed_profile_tx = [];
%     for k=1:tx
%         speed_profile_tx(:,:,k) = fftshift(fft(range_profile_tx(:,:,k), M, 2));
%     end
%     
%     figure(4)
%     imagesc(abs(squeeze(speed_profile(:,:,1).')))
%     figure(5)
%     imagesc(abs(squeeze(fftshift(range_profile(:,:,1).', 1))))
%     
%     figure(1)
%     subplot(221)
%     mesh(abs(squeeze(speed_profile(:,:,1))))
%     subplot(222)
%     mesh(abs(squeeze(speed_profile_tx(:,:,1))))
%     
% %     CACFAR
%     detect_profile = [];
%     detect_res = {}; % range, doppler, val
%     for k=1:n_RX
%         [detect_profile(:,:,k), detect_res{k}] = CFAR(squeeze(speed_profile(:,:,k)));
%     end
%     
%     detect_profile_tx = [];
%     detect_res_tx = {}; % range, doppler, val
%     for k=1:tx
%         [detect_profile_tx(:,:,k), detect_res_tx{k}] = CFAR(squeeze(speed_profile_tx(:,:,k)));
%     end
% 
%     detect_filtered = DetectPointFilter(detect_res);
%     detect_filtered_tx = DetectPointFilter(detect_res_tx);
%     
%     angle_az = [];
%     angle_el = [];
%     range = [];
%     %% 方位向角度
%     for k=1:size(detect_filtered,2)
%         ii = detect_filtered(1, k);
%         jj = detect_filtered(2, k);
%         
%         angle_data = squeeze(speed_profile(ii,jj,:));
%         angle_fft = fft(angle_data, Q);
%         [iq,index] = max(abs(angle_fft));
%         angle_az = [angle_az an_axis_az(index)];
%         
%         % 俯仰角度
%         angle_data = squeeze(speed_profile_tx(ii,jj,:));
%         angle_fft = fft(angle_data, Q);
%         [iq,index] = max(abs(angle_fft));
%         angle_el = [angle_el an_axis_el(index)];
%         
%         range = [range range_axis(ii)];
%     end
% 
%     subplot(223)
%     scatter(range .* cos(angle_az), range .* sin(angle_az), 'x')
%     axis([-10 10 0 5])
%     subplot(224)
%     scatter(range .* cos(angle_el), range .* sin(angle_el), 'x')
%     axis([-10 10 0 5])
%     
%     figure(2)
%     range_xy = range .* cos(angle_el);
%     range_z = range .* sin(angle_el);
%     range_x = range_xy .* cos(angle_az);
%     range_y = range_xy .* sin(angle_az);
%     scatter3(range_x, range_y, range_z);
%     axis([0, 10, -5, 5, -5, 5])
%     view(-88,10)
%     az_rmap = fftshift(squeeze(sum(abs(fft(speed_profile, Q, 3)), 2)), 2);
%     el_rmap = fftshift(squeeze(sum(abs(fft(speed_profile_tx, Q, 3)), 2)), 2);
%     
%     angle_detect_az = AnglePeaksDetect(az_rmap);
%     angle_detect_el = AnglePeaksDetect(el_rmap);
%     
%     if isempty(angle_detect_az) || isempty(angle_detect_el)
%         continue;
%     end
%     
%     subplot(223)
%     scatter(range_axis(angle_detect_az(:, 1)) .* cos(an_axis_az(angle_detect_az(:, 2))), ...
%         range_axis(angle_detect_az(:, 1)) .* sin(an_axis_az(angle_detect_az(:, 2))), 'x')
%     axis([-20 20 0 25])
%     subplot(224)
%     scatter(range_axis(angle_detect_el(:, 1)) .* cos(an_axis_el(angle_detect_el(:, 2))), ...
%         range_axis(angle_detect_el(:, 1)) .* sin(an_axis_el(angle_detect_el(:, 2))), 'x')
%     axis([-20 20 0 25])
%     
%     extract_angle_range = AngleSearchEqual(angle_detect_az, angle_detect_el);
%     
%     if ~isempty(extract_angle_range)
%         range = range_axis(extract_angle_range(:,1));
%         angle_az = an_axis_az(extract_angle_range(:,2));
%         angle_el = an_axis_el(extract_angle_range(:,3));
% 
%         figure(2)
%         range_xy = range .* cos(angle_el);
%         range_z = range .* sin(angle_el);
%         range_x = range_xy .* cos(angle_az);
%         range_y = range_xy .* sin(angle_az);
%         scatter3(range_x, range_y, range_z);
%         axis([0, 25, -20, 20, -15, 15])
%         view(-80,30)
%     end
%     
%     figure(3)
%     mesh(az_rmap);
%     view(-1,80)
%     
%     figure(4)
%     mesh(el_rmap);
%     view(-1,80)
%     
%     drawnow

end

fclose(fid);
