clear all;
close all;
clc;

addpath algorithm
addpath building
addpath util

wallStart = [-1.0, 2.2];
wallEnd   = [3.0, 2.2];
barrierStart = [2.2 1.7];
barrierEnd   = [2.2 -1.0];
humanStart = [4.2, 1];
humanEnd   = [4.2, -1];

GHz = 1e9;
MHz = 1e6;
us = 1e-6;
%% 雷达参数（使用mmWave Studio默认参数）
c=3.0e8;  
B=4000*1e6;       %调频带宽
K = 68*1e12;      %调频斜率
T=B/K;            %采样时间
Tc=65e-6;         %chirp总周期
fs=10*1e6;       %采样率
f0=77e9;          %初始频率
lambda=c/f0;      %雷达信号波长
d=lambda/2;       %天线阵列间距
n_samples=256;    %采样点数/脉冲
N=256;            %距离向FFT点数
n_chirps=64;      %每帧脉冲数
M=128;            %多普勒向FFT点数
n_RX=4;           %RX天线通道数
n_TX=2;           %TX天线通道数
Q = 512;          %角度FFT
tx = 2;           %发射天线数目
rx = 4;          %接受天线数目
fnumber = 256;
PRI = 4e-3;

emitterAngle  = pi / 4;
correctMatrix = [cos(emitterAngle) -sin(emitterAngle); sin(emitterAngle) cos(emitterAngle)];
fontSize = 16;

range_axis = linspace(0, fs * c / 2 / K, N);
velo_axis = linspace(-lambda / 4 / Tc, lambda / 4 / Tc, M);
an_axis_az = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q); % 雷达是倒立过来的
an_axis_el = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
fname = "2023-05-25-15-34-03.bin";
xxx = [];
yyy = [];

% 探测结果
confirm_detect = [];
temp_detect = [];
origin_detect = []; % 不加入任何抗干扰措施的结果


fid = fopen(fname,'rb'); 
sdata = fread(fid,'int16');  
tx_select_rx = 1;
rx_select_tx = 1;  
    tic
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
    
    
%     t_lvds_data = [];
    for tr = 1:tx * rx
        t_lvds_data(tr, :, :) = lvds_data((tr - 1) * fnumber+1:tr * fnumber, :);
    end
    
%     tx_data = [];
    
%     data_tx = [];
    
    data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
    data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
    data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
    data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));       %RX4

%     for ttxx = 1:tx
%         if ttxx == rx_select_tx
%             data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
%             data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
%             data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
%             data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));       %RX4
%         end
%         
% %         temp_data = squeeze(t_lvds_data(ttxx * rx - 3, :, :));          % 采用RX1
% %         data_tx(:, :, ttxx) = temp_data(:,1:end-8) - temp_data(:,9:end);
%     end
   
    data_radar_11 = data_radar_1(:,1:end-8) - data_radar_1(:,9:end);
    data_radar_22 = data_radar_2(:,1:end-8) - data_radar_2(:,9:end);
    data_radar_33 = data_radar_3(:,1:end-8) - data_radar_3(:,9:end);
    data_radar_44 = data_radar_4(:,1:end-8) - data_radar_4(:,9:end);
    
%     data_radar=[];            
    data_radar(:,:,1) = data_radar_11;     %三维雷达回波数据
    data_radar(:,:,2) = data_radar_22;
    data_radar(:,:,3) = data_radar_33;
    data_radar(:,:,4) = data_radar_44;
    %% 3维FFT处理
    %距离FFT
%     range_win = hamming(n_samples);   %加海明窗
%     doppler_win = hamming(n_chirps);
%     range_profile = [];
    for k=1:n_RX
        range_profile(:,:,k) = fft(data_radar(:,:,k), N);
%         speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    speed_profile = fftshift(fft(range_profile(:,:,1), M, 2), 2);
    az_rmap = fftshift(squeeze(sum(abs(fft(range_profile, Q, 3)), 2)), 2);    % 距离 x 角度维度
    
%     range_profile_tx = [];
%     for k=1:tx
%         range_profile_tx(:,:,k) = fft(data_tx(:,:,k), N);
%     end
    
    %多普勒FFT
%     speed_profile = [];
%     for k=1:n_RX
%         speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
%     end
   
%     speed_profile_tx = [];
%     for k=1:n_TX
%         speed_profile_tx(:,:,k) = fftshift(fft(range_profile_tx(:,:,k), M, 2), 2);
%     end
    
%     rd_result = squeeze(sum(speed_profile, 3) / size(speed_profile, 3));

    [detect_profile, detect_res] = CFAR(speed_profile);
    %% detect_res(:, 2) is range 
%     el_rmap = fftshift(squeeze(sum(abs(fft(speed_profile_tx, Q, 3)), 2)), 2); % 距离 x 角度维度
    angle_detect = AngleDetectSearchingAngle(az_rmap, detect_res);
    % 请用 (detect_res(:, 2), detect_res(:, 1) + 9) 来标定 rd 谱图与检测 坐标轴
    % detect_res(:, 1) 是 距离轴
    if isempty(angle_detect)
        xxxx = [];
        yyyy = [];
    else
        angles = an_axis_az(angle_detect(:,2));
        ranges = range_axis(angle_detect(:,1));
%         detect_result_UOI = [ranges; angles];
%         
%         drawUOIarea(az_rmap, detect_result_UOI, range_axis, an_axis_az);

        xxxx = ranges .* sin(angles);
        yyyy = ranges .* cos(angles);
    end
    
    origin_detect = [origin_detect; xxxx' yyyy'];
    % 
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
%         if ~isempty(new_detect_res)
%             xxx = [xxx new_detect_res(:, 1).'];
%             yyy = [yyy new_detect_res(:, 2).'];
%         end

% figure(10006)
%     suptitle("Radar waveform（Real Part）");
%     subplot(221)
%     plot(real(data_radar_1(:, 1)))
%     axis([0 256 min(real(data_radar_1(:, 1))) max(real(data_radar_1(:, 1)))])
%     subplot(222)
%     plot(real(data_radar_2(:, 1)))
%     axis([0 256 min(real(data_radar_2(:, 1))) max(real(data_radar_2(:, 1)))])
%     subplot(223)
%     plot(real(data_radar_3(:, 1)))
%     axis([0 256 min(real(data_radar_3(:, 1))) max(real(data_radar_3(:, 1)))])
%     subplot(224)
%     plot(real(data_radar_4(:, 1)))
%     axis([0 256 min(real(data_radar_4(:, 1))) max(real(data_radar_4(:, 1)))])
% figure(10007)
% %     range_profile1 = flip(fft(data_radar_1, N));
%     imagesc(linspace(1, 64), range_axis, abs(squeeze(range_profile(:,:,1))))
%     xlabel("Chirps", 'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance (m)", 'fontsize',fontSize,'FontName','Times New Roman')
%     title("Range Profile", 'fontsize',fontSize,'FontName','Times New Roman')
% figure(10008)
% %     speed_profile1 = fft(range_profile1, M, 2);
%     imagesc(velo_axis, range_axis, abs(speed_profile(:,:,1)))
%     xlabel("velocity (m/s)", 'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance (m)", 'fontsize',fontSize,'FontName','Times New Roman')
%     title("Range-Doppler Result", 'fontsize',fontSize,'FontName','Times New Roman')
% figure(10009)
%     imagesc(velo_axis, range_axis, abs(detect_profile))
%     xlabel("velocity (m/s)", 'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance (m)", 'fontsize',fontSize,'FontName','Times New Roman')
%     title("CA-CFAR Result", 'fontsize',fontSize,'FontName','Times New Roman')
% figure(10010)
%     imagesc(an_axis_az, range_axis, abs(az_rmap));
%     xlabel("degree(°)",'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance(m)",'fontsize',fontSize,'FontName','Times New Roman')
%     title("Azimuth Range-Angle Result",'fontsize',fontSize,'FontName','Times New Roman')
% figure(10011)
%     scatter(xxx, yyy, 'filled')
%     axis([-6 6 0 10])
%     xlabel("X axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("Y axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
%     title("Locate results",'fontsize',fontSize,'FontName','Times New Roman')
%     hold on
%     plot_human_trace([4.6 3.2], [4.6 5.2], pi / 4)
%     plot_wall_trace(wallStart, wallEnd, pi / 4)
%     plot_wall_trace(barrierStart, barrierEnd, pi / 4)
%     legend('Locate results', 'Actual trajectory', 'Steel wall', 'Wooden board')
% figure(10012)
%     imagesc(an_axis_el, range_axis, abs(el_rmap));
%     xlabel("degree(°)",'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance(m)",'fontsize',fontSize,'FontName','Times New Roman')
%     title("Elevation Range-Angle Result",'fontsize',fontSize,'FontName','Times New Roman')
% figure(10013)
%     imagesc(db(az_rmap(1:55,:)));
%     xlabel("degree(°)",'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance(m)",'fontsize',fontSize,'FontName','Times New Roman')
%     title("Azimuth Range-Angle Result",'fontsize',fontSize,'FontName','Times New Roman')
% figure(10014)
%     imagesc(db(speed_profile(1:55,:,1)))
%     xlabel("velocity (m/s)", 'fontsize',fontSize,'FontName','Times New Roman')
%     ylabel("distance (m)", 'fontsize',fontSize,'FontName','Times New Roman')
%     title("Range-Doppler Result", 'fontsize',fontSize,'FontName','Times New Roman')
end    
toc

cc = [xxx' yyy'];
if ~isempty(cc)
    figure(10011)
    scatter(xxx, yyy, 'filled')
    axis([-6 6 0 10])
    xlabel("X axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
    ylabel("Y axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
    title("Locate results",'fontsize',fontSize,'FontName','Times New Roman')
    hold on
    plot_human_trace([4.6 3.2], [4.6 5.2], pi / 4)
    plot_wall_trace(wallStart, wallEnd, pi / 4)
    plot_wall_trace(barrierStart, barrierEnd, pi / 4)
    legend('Locate results', 'Actual trajectory', 'Steel wall', 'Wooden board')

    
    [correct_points, mirror_loc] = wallMultipathCorrect(pi / 4, wallStart, wallEnd, cc);

    figure(11000)
    scatter(correct_points(:,1), correct_points(:,2),'filled')
    axis([-6 6 -1.5 5])
    hold on
    drawWall(wallStart, wallEnd, 3)
    drawWall(barrierStart, barrierEnd, 2)
    drawPersonTrace([4.2, 1], [4.2, -1])
    legend('Detect result', 'Wall', 'Board', 'Real trajectory')
    xlabel('X axis (m)','fontsize',fontSize,'FontName','Times New Roman');
    ylabel('Y axis (m)','fontsize',fontSize,'FontName','Times New Roman');
    title("The correct result",'fontsize',fontSize,'FontName','Times New Roman');
else
    fprintf("无探测结果。");
end

fclose(fid);
