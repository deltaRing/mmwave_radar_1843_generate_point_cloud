clear all;
close all;
clc;

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
B=1800*1e6;       %调频带宽
K = 68*1e12;      %调频斜率
T=B/K;            %采样时间
Tc=65e-6;         %chirp总周期
fs=5*1e6;         %采样率
f0=77e9;          %初始频率
lambda=c/f0;      %雷达信号波长
d=lambda/2;       %天线阵列间距
n_samples=256;    %采样点数/脉冲
N=256;            %距离向FFT点数
n_chirps=128;     %每帧脉冲数
M=128;            %多普勒向FFT点数
n_RX=4;           %RX天线通道数
n_TX=3;           %TX天线通道数
Q = 512;          %角度FFT
tx = 3;           %发射天线数目
rx = 4;          %接受天线数目
rmax = fs * c / 2 / K;
fnumber = 256;
PRI = 4e-3;

emitterAngle  = pi / 4;
correctMatrix = [cos(emitterAngle) -sin(emitterAngle); sin(emitterAngle) cos(emitterAngle)];
fontSize = 16;

range_axis = linspace(0, fs * c / 2 / K, N);
velo_axis = linspace(-lambda / 4 / Tc, lambda / 4 / Tc, M);
an_axis_az = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q); % 雷达是倒立过来的
an_axis_el = linspace(-asin(lambda/2/d), asin(lambda/2/d), Q);
fname = "T3R4_f500_rate5000.bin";
xxx = [];
yyy = [];

% 探测结果
confirm_detect = [];
temp_detect = [];
origin_detect = []; % 不加入任何抗干扰措施的结果

antenna_loc = [lambda * 3 / 2 0; lambda 0; lambda / 2 0; 0 0];
range_max = fs * c / 2 / K;

fid = fopen(fname,'rb'); 
sdata = fread(fid,'int16');  
tx_select_rx = 1;
rx_select_tx = 1;    
    tic
for xx = 1:fnumber-1 
    tic
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
    
    for tr = 1:tx * rx
        t_lvds_data(tr, :, :) = lvds_data((tr - 1) * fnumber+1:tr * fnumber, :);
    end
    
    range_win = hamming(n_samples);   %加海明窗
    data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :)); % .* range_win;   %RX1
    data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :)); % .* range_win;   %RX2
    data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :)); % .* range_win;   %RX3
    data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :)); % .* range_win;       %RX4
    % data_radar_11 = data_radar_1(:,1:end-1) - data_radar_1(:,2:end);
    % data_radar_22 = data_radar_2(:,1:end-1) - data_radar_2(:,2:end);
    % data_radar_33 = data_radar_3(:,1:end-1) - data_radar_3(:,2:end);
    % data_radar_44 = data_radar_4(:,1:end-1) - data_radar_4(:,2:end);
    data_radar_11 = data_radar_1;
    data_radar_22 = data_radar_2;
    data_radar_33 = data_radar_3;
    data_radar_44 = data_radar_4;
    
%     data_radar=[];            
    data_radar(:,:,1) = data_radar_11;     %三维雷达回波数据
    data_radar(:,:,2) = data_radar_22;
    data_radar(:,:,3) = data_radar_33;
    data_radar(:,:,4) = data_radar_44;
    %% 3维FFT处理
    %距离FFT
%     doppler_win = hamming(n_chirps);
%     range_profile = [];
    for k=1:n_RX
        range_profile(:,:,k) = fft(data_radar(:,:,k), N);
%         speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    speed_profile = fftshift(fft(range_profile(:,:,1), M, 2), 2);
    az_rmap = fftshift(squeeze(sum(abs(fft(range_profile, Q, 3)), 2)), 2);    % 距离 x 角度维度
    
    data = squeeze(range_profile(:,1,:));
    
    [detect_profile, detect_res] = CA_CFAR(speed_profile);
    res = RangeCentroid(detect_res);
    %% detect_res(:, 2) is range 
    [range, angle] = getAngleInfo(res, az_rmap);
    % 请用 (detect_res(:, 2), detect_res(:, 1) + 9) 来标定 rd 谱图与检测 坐标轴
    % detect_res(:, 1) 是 距离轴
    if isempty(range)
        xxxx = [];
        yyyy = [];
    else
        xxxx = [];
        yyyy = [];
        for aaa = 1:length(angle) 
            angles = an_axis_az(angle{aaa});
            ranges = range_axis(range(aaa));

            xxxx = [xxxx ranges .* sin(angles)];
            yyyy = [yyyy ranges .* cos(angles)];
        end
    end
    % 不滤波的
    xxx = [xxx xxxx];
    yyy = [yyy yyyy];

    % 滤波的
    % origin_detect = [origin_detect; xxxx' yyyy'];
    % % 
    % new_detect_res = zeros(length(xxxx), 3);
    % new_detect_res(:, 1) = xxxx; 
    % new_detect_res(:, 2) = yyyy;
    % % 填充数据
    % 
    % [detect_res_, temp_detect_res] = PointTraceFilter(xx,...
    %     new_detect_res, ...
    %     confirm_detect, ...
    %     temp_detect);
    % 
    % confirm_detect = detect_res_;
    % temp_detect = temp_detect_res;
    % 
    % if ~isempty(confirm_detect)
    %     xxx = [xxx confirm_detect(:, 1).'];
    %     yyy = [yyy confirm_detect(:, 2).'];
    % end

    % 
    % figure(12001)
    % imagesc(velo_axis, range_axis, abs(speed_profile))
    % xlabel("velocity (m/s)", 'fontsize',fontSize,'FontName','Times New Roman')
    % ylabel("distance (m)", 'fontsize',fontSize,'FontName','Times New Roman')
    % title("Range-Doppler Result", 'fontsize',fontSize,'FontName','Times New Roman')

    figure(12000)
    scatter(xxx, yyy, 5, 'filled')
    axis([-12.5 12.5 -1.5 15])
    drawnow
    % 填充数据
    toc
end    
toc


fclose(fid);