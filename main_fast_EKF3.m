% main_fast _EKF2_
% 试图加入EKF来增加NLOS目标行走的稳定性
% main_fast_EKF3.m
% 不要动任何参数！谁知道会发生什么事

clear all;
close all;
clc;

addpath trace_optimize
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
fname = "2023-05-25-15-26-25.bin";

% 记录所提方法的定位结果
xxx = []; yyy = []; zzz = [];
xxx2 = []; yyy2 = []; zzz2 = [];

% 探测结果
confirm_detect = [];
temp_detect = [];
origin_detect = []; % 不加入任何抗干扰措施的结果

% 启用俯仰角信息
ELevation_Enable = 1;
% 记录trace_id
Trace_ID = 0;

fid = fopen(fname,'rb'); 
sdata = fread(fid,'int16');  
tx_select_rx = 1;
rx_select_tx = 1;   
el_channel   = 2;
% 选择的俯仰向通道

% 初始化协方差 Pn
% 和运动状态 Xn
Pn = {};
Xn = [];

Xnn = [];

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
    
    for tr = 1:tx * rx
        t_lvds_data(tr, :, :) = lvds_data((tr - 1) * fnumber+1:tr * fnumber, :);
    end
    
    data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
    data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
    data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
    data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));       %RX4
    data_radar_5 = squeeze(t_lvds_data(el_channel   * rx - 3, :, :));   %RX1_TX2
    data_radar_11 = data_radar_1(:,1:end-8) - data_radar_1(:,9:end);
    data_radar_22 = data_radar_2(:,1:end-8) - data_radar_2(:,9:end);
    data_radar_33 = data_radar_3(:,1:end-8) - data_radar_3(:,9:end);
    data_radar_44 = data_radar_4(:,1:end-8) - data_radar_4(:,9:end);
    data_radar_55 = data_radar_5(:,1:end-8) - data_radar_5(:,9:end);
    
    
%     data_radar=[];            
    data_radar(:,:,1) = data_radar_11;     %三维雷达回波数据
    data_radar(:,:,2) = data_radar_22;
    data_radar(:,:,3) = data_radar_33;
    data_radar(:,:,4) = data_radar_44;
    data_radar(:,:,5) = data_radar_55;
    %% 3维FFT处理
    %距离FFT
%     range_win = hamming(n_samples);   %加海明窗
%     doppler_win = hamming(n_chirps);
%     range_profile = [];
    for k=1:n_RX + 1
        range_profile(:,:,k) = fft(data_radar(:,:,k), N);
        speed_profile(:,:,k) = fftshift(fft(range_profile(:,:,k), M, 2), 2);
    end
    range_profile_az = [range_profile(:,:,1:n_RX)]; % 方位向的距离像
    range_profile_el = [range_profile(:,:,[1 5])];   % 俯仰向的距离像
    speed_profile1 = fftshift(fft(range_profile(:,:,1), M, 2), 2);
    az_rmap = fftshift(squeeze(sum(abs(fft(range_profile_az, Q, 3)), 2)), 2);    % 距离 x 角度维度
    el_rmap = fftshift(squeeze(sum(abs(fft(range_profile_el, Q, 3)), 2)), 2);    % 距离 x 角度维度

    figure(1)
    imagesc(db(el_rmap))

    [detect_profile, detect_res] = CFAR(speed_profile1);
    %% detect_res(:, 2) is range 
    angle_detect = AngleDetectSearchingAngle(az_rmap, detect_res);
    % 请用 (detect_res(:, 2), detect_res(:, 1) + 9) 来标定 rd 谱图与检测 坐标轴
    % detect_res(:, 1) 是 距离轴
    if isempty(angle_detect)
        xxxx = [];
        yyyy = [];
        zzzz = [];
    else
        angles = an_axis_az(angle_detect(:,2));
        if ELevation_Enable % 俯仰角
            angle_elevation      = el_rmap(angle_detect(:,1), :);
            [~, elevation_index] = max(abs(angle_elevation).');
            elevation            = an_axis_el(elevation_index); 
            ranges               = range_axis(angle_detect(:,1));
            ranges_xy            = ranges .* cos(elevation);
            xxxx                 = ranges_xy .* sin(angles);
            yyyy                 = ranges_xy .* cos(angles);
            zzzz                 = ranges .* sin(elevation);
            xxx2 = [xxx2 xxxx]; yyy2 = [yyy2 yyyy]; zzz2 = [zzz2 zzzz];
        else
            ranges = range_axis(angle_detect(:,1));
            xxxx = ranges .* sin(angles);
            yyyy = ranges .* cos(angles);
            xxx2 = [xxx2 xxxx]; yyy2 = [yyy2 yyyy];
        end

    end
    
    if ELevation_Enable
        origin_detect = [origin_detect; xxxx' yyyy' zzzz'];
        new_detect_res = zeros(length(xxxx), 3);
        new_detect_res(:, 1) = xxxx; 
        new_detect_res(:, 2) = yyyy;
        new_detect_res(:, 3) = zzzz;
    else
        origin_detect = [origin_detect; xxxx' yyyy'];
        new_detect_res = zeros(length(xxxx), 3);
        new_detect_res(:, 1) = xxxx; 
        new_detect_res(:, 2) = yyyy;
    end

    % 填充数据
    [detect_res_, temp_detect_res, Trace_ID] = PointTraceFilter(xx,...
        new_detect_res, ...
        confirm_detect, ...
        temp_detect, ...
        Trace_ID);
        
    % 检查是否有新的航迹
    [res, new_index, remove_index, loss_index] = check_new_trace(detect_res_, confirm_detect);
    % 先删除
    if ~isempty(remove_index)
        [Xn, Pn] = remove_EKF(Xn, Pn, remove_index);
    end
    % 初始化
    if ~isempty(new_index)
        [Xn, Pn] = init_EKF(detect_res_(new_index, 1), ...
                    detect_res_(new_index, 2), ...
                    detect_res_(new_index, 3), ...
                    Xn, Pn);
    end
    % 再来EKF
    if ~isempty(Pn)
        Aindex = 1;
        A = {}; E = {}; R = {};
        for ii = 1:size(detect_res_, 1)
            if ~isempty(find(loss_index == ii))
                A{Aindex} = {};
                E{Aindex} = {};
                R{Aindex} = {};
                Aindex = Aindex + 1;
            else
                px = detect_res_(ii, 1); py = detect_res_(ii, 2); pz = detect_res_(ii, 3);
                a = atan2(py, px);
                e = atan2(pz, norm([px py]));
                r = norm([px py pz]);
                A{Aindex} = a;
                E{Aindex} = e;
                R{Aindex} = r;
                Aindex = Aindex + 1;
            end
        end

        [Xn, Pn] = EKF3(E, A, R, Xn, Pn);
        for ii = 1:size(Xn,1)
            detect_res_(ii, 1) = Xn(ii, 1);
            detect_res_(ii, 2) = Xn(ii, 3);
            detect_res_(ii, 3) = Xn(ii, 5);
        end
    end

    confirm_detect = detect_res_;
    temp_detect = temp_detect_res;
    
    if ~isempty(confirm_detect)
        if ELevation_Enable 
            xxx = [xxx confirm_detect(:, 1).'];
            yyy = [yyy confirm_detect(:, 2).'];
            zzz = [zzz confirm_detect(:, 3).'];
        else
            xxx = [xxx confirm_detect(:, 1).'];
            yyy = [yyy confirm_detect(:, 2).'];
        end
    end
end    
toc

cc = [xxx' yyy'];
if ~isempty(cc)
    figure(10011)
    scatter3(xxx, yyy, zzz, 'filled')
    % axis([-6 6 0 10])
    xlabel("X axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
    ylabel("Y axis (m)",'fontsize',fontSize,'FontName','Times New Roman')
    title("Locate results using the proposed method",'fontsize',fontSize,'FontName','Times New Roman')
    hold on
    plot_human_trace([4.6 3.2], [4.6 5.2], pi / 4)
    plot_wall_trace(wallStart, wallEnd, pi / 4)
    plot_wall_trace(barrierStart, barrierEnd, pi / 4)
    legend('Locate results', 'Actual trajectory', 'Steel wall', 'Wooden board')

    % [correct_points, mirror_loc] = wallMultipathCorrect(pi / 4, wallStart, wallEnd, cc);
    % figure(11000)
    % scatter3(correct_points(:,1), correct_points(:,2), zzz,'filled')
    % axis([-6 6 -1.5 5])
    % hold on
    % drawWall(wallStart, wallEnd, 3)
    % drawWall(barrierStart, barrierEnd, 2)
    % drawPersonTrace([4.2, 1], [4.2, -1])
    % legend('Detect result', 'Wall', 'Board', 'Real trajectory')
    % xlabel('X axis (m)','fontsize',fontSize,'FontName','Times New Roman');
    % ylabel('Y axis (m)','fontsize',fontSize,'FontName','Times New Roman');
    % title("The correct result",'fontsize',fontSize,'FontName','Times New Roman');
else
    fprintf("无探测结果。");
end

fclose(fid);