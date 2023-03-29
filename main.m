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
M=512;         %多普勒向FFT点数
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
       lvds_data(1,count) =1i* sdata2(i) + sdata2(i+2); 
       lvds_data(1,count+1) =1i* sdata2(i+1)+sdata2(i+3); 
       count = count + 2;
    end
    lvds_data = reshape(lvds_data, tx*n_samples*n_RX, n_chirps);
    
    
    t_lvds_data = [];
    for tr = 1:tx * rx
        t_lvds_data(tr, :, :) = lvds_data((tr - 1) * fnumber+1:tr * fnumber, :);
    end
    
    tx_data = [];
    tx_select_rx = 1;
    rx_select_tx = 1;
    for ttxx = 1:tx
%         lvds_data_tx = squeeze(lvds_data(ttxx,:,:));
%         lvds_data_tx = lvds_data_tx.';
%         cdata = zeros(n_RX,n_chirps*n_samples);
%         for row = 1:n_RX
%               for i = 1: n_chirps
%                   cdata(row,(i-1)*n_samples+1:i*n_samples) = lvds_data_tx(i,(row-1)*n_samples+1:row*n_samples);
%               end
%               if tx_select_rx == row
%                   temp_data = reshape(cdata(row,:),n_samples,n_chirps);
% %                   temp_data = temp_data(:,2:end) - temp_data(:,1:end-1);
%                   tx_data(:,:,ttxx) = temp_data;
%               end
%         end
%         
%         if ttxx == rx_select_tx
%             data_radar_1 = reshape(cdata(1,:),n_samples,n_chirps);   %RX1
%             data_radar_2 = reshape(cdata(2,:),n_samples,n_chirps);   %RX2
%             data_radar_3 = reshape(cdata(3,:),n_samples,n_chirps);   %RX3
%             data_radar_4 = reshape(cdata(4,:),n_samples,n_chirps);   %RX4
%         end

        if ttxx == rx_select_tx
            data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
            data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
            data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
            data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));   %RX4
        end
    end
    
    figure(6)
    subplot(221)
    plot(real(squeeze(data_radar_1(:,1))))
    subplot(222)
    plot(real(squeeze(data_radar_2(:,1))))
    subplot(223)
    plot(real(squeeze(data_radar_3(:,1))))
    subplot(224)
    plot(real(squeeze(data_radar_4(:,1))))
    
    data_radar_1 = data_radar_1(:,7:end) - data_radar_1(:,1:end-6);
    data_radar_2 = data_radar_2(:,7:end) - data_radar_2(:,1:end-6);
    data_radar_3 = data_radar_3(:,7:end) - data_radar_3(:,1:end-6);
    data_radar_4 = data_radar_4(:,7:end) - data_radar_4(:,1:end-6);
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
        range_profile(:,:,k) = fft(data_radar(:,:,k), N);
    end
    %多普勒FFT
    speed_profile = [];
    for k=1:n_RX
        speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    
%     range_profile1 = flip(fft(data_radar_1, N));
%     figure(10001)
%     imagesc(abs(squeeze(range_profile(:,:,1))))
% %     speed_profile1 = fft(range_profile1, M, 2);
%     figure(10002)
%     imagesc(abs(speed_profile(:,:,1)))
    
 %     CACFAR
    [detect_profile, detect_res] = CFAR(squeeze(speed_profile(:,:,1)));
    %% detect_res(:, 2) is range 
    
%     figure(10003)
%     mesh(abs(squeeze(speed_profile(:,:,1))))
%     figure(10004)
%     mesh(abs(detect_profile))


    az_rmap = squeeze(sum(abs(fft(speed_profile, Q, 3)), 2));
    angle_detect = AnglePeaksDetect(az_rmap);
%     figure(10005)
%     mesh(abs(az_rmap));
    
    angles = an_axis_az(angle_detect(:,2));
    ranges = range_axis(angle_detect(:,1));
    
    xx = ranges .* sin(angles);
    yy = ranges .* cos(angles);
    
    figure(10006)
    scatter(xx, yy)
    hold on
    axis([-15 15 0 20])

figure(10007)
subplot(231)
    title("距离像结果")
    range_profile1 = flip(fft(data_radar_1, N));
    imagesc(abs(squeeze(range_profile(:,:,1))))
subplot(232)
    title("RD谱图结果")
    speed_profile1 = fft(range_profile1, M, 2);
    imagesc(abs(speed_profile(:,:,1)))
subplot(233)
    title("CA-CFAR结果")
    mesh(abs(detect_profile))
subplot(234)
    title("方位向距离-角度谱图结果")
    mesh(abs(az_rmap));
subplot(235)
    title("定位结果")
    scatter(xx, yy)
    hold on
    axis([-15 15 0 20])
subplot(236)

end

    fclose(fid);
