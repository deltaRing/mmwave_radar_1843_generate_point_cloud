clear all;close all;clc;

addpath algorithm

GHz = 1e9;
MHz = 1e6;
us = 1e-6;
%% 雷达参数（使用mmWave Studio默认参数）
c=3.0e8;  
B=4000*1e6;       %调频带宽
K = 58*1e12;  %调频斜率
T=B/K;         %采样时间
Tc=65e-6;     %chirp总周期
fs=5.5*1e6;       %采样率
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
PRI = 4e-3;

range_axis = linspace(0, fs * c / 2 / K, N);
velo_axis = linspace(-lambda / 4 / Tc, lambda / 4 / Tc, M);
an_axis_az = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
an_axis_el = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
fname = "2023-05-25-15-34-03.bin";
xxx = [];
yyy = [];

% 探测结果
confirm_detect = [];
temp_detect = [];


fid = fopen(fname,'rb'); 
sdata = fread(fid,'int16');    
for xx = 1:fnumber-1   
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
    
    data_tx = [];
    
    for ttxx = 1:tx
        if ttxx == rx_select_tx
            data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
            data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
            data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
            data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));   %RX4
        end
        
        data_tx(:, :, ttxx) = data_radar_1(:,7:end) - data_radar_1(:,1:end-6);
    end
   
    data_radar_1 = data_radar_1(:,1:end-8) - data_radar_1(:,9:end);
    data_radar_2 = data_radar_2(:,1:end-8) - data_radar_2(:,9:end);
    data_radar_3 = data_radar_3(:,1:end-8) - data_radar_3(:,9:end);
    data_radar_4 = data_radar_4(:,1:end-8) - data_radar_4(:,9:end);
    
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
    
    range_profile_tx = [];
    for k=1:tx
        range_profile_tx(:,:,k) = fft(data_tx(:,:,k), N);
    end
    
    %多普勒FFT
    speed_profile = [];
    for k=1:n_RX
        speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    
    rd_result = squeeze(sum(speed_profile, 3) / size(speed_profile, 3));
    [detect_profile, detect_res] = CFAR(rd_result);
    %% detect_res(:, 2) is range 
    az_rmap = fftshift(squeeze(sum(abs(fft(speed_profile, Q, 3)), 2)), 2); % 距离 x 角度维度
%     el_rmap = squeeze(sum(abs(fft(speed_profile_tx, Q, 3)), 2));
    angle_detect = AngleDetectSearchingAngle(az_rmap, detect_res);
    % 请用 (detect_res(:, 2), detect_res(:, 1) + 9) 来标定 rd 谱图与检测 坐标轴
    % detect_res(:, 1) 是 距离轴
    if isempty(angle_detect)
        xxxx = [];
        yyyy = [];
    else
        angles = an_axis_az(angle_detect(:,2));
        ranges = range_axis(angle_detect(:,1));

        xxxx = ranges .* sin(angles);
        yyyy = ranges .* cos(angles);
    end
    
    new_detect_res = zeros(length(xxxx), 3);
    new_detect_res(:, 1) = xxxx; 
    new_detect_res(:, 2) = yyyy;
    % 填充数据
    
    [detect_res_, temp_detect_res] = PointTraceFilter(xx,...
        new_detect_res, ...
        confirm_detect, ...
        temp_detect);
        
    confirm_detect = detect_res_;
    temp_detect = temp_detect_res;
    
    if ~isempty(confirm_detect)
        xxx = [xxx confirm_detect(:, 1).'];
        yyy = [yyy confirm_detect(:, 2).'];
    end

figure(10007)
    range_profile1 = flip(fft(data_radar_1, N));
    imagesc(linspace(1, 64), range_axis, abs(squeeze(range_profile(:,:,1))))
    xlabel("Chirps")
    ylabel("距离(m)")
    title("距离像结果")
figure(10008)
    speed_profile1 = fft(range_profile1, M, 2);
    imagesc(velo_axis, range_axis, abs(speed_profile(:,:,1)))
    xlabel("速度(m/s)")
    ylabel("距离(m)")
    title("RD谱图结果")
figure(10009)
    imagesc(velo_axis, range_axis, abs(detect_profile))
    xlabel("速度(m/s)")
    ylabel("距离(m)")
    title("CA-CFAR结果")
figure(10010)
    imagesc(an_axis_az, range_axis, abs(az_rmap));
    xlabel("度(°)")
    ylabel("距离(m)")
    title("方位向距离-角度谱图结果")
figure(10011)
    scatter(xxx, yyy)
    axis([-15 15 0 20])
    xlabel("X轴 距离(m)")
    ylabel("Y轴 距离(m)")
    title("定位结果")
figure(10012)
    mesh(abs(az_rmap));
    xlabel("度(°)")
    ylabel("距离(m)")
    title("方位向距离-角度谱图结果")

end

    fclose(fid);
