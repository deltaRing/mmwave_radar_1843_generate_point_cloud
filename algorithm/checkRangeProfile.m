% 查看距离像的脚本
function checkRangeProfile(filename, fnumber)
    if nargin == 1
        fnumber = 256;
    end
    if nargin == 0
        error("至少需要参数：文件名字，可选参数：fnumber 帧数");
        return;
    end

    fid = fopen(filename,'rb'); 
    sdata = fread(fid,'int16'); 
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
    M=256;         %多普勒向FFT点数
    n_RX=4;        %RX天线通道数
    Q = 512;       %角度FFT
    tx = 2;        %发射天线数目
    rx = 4;        %接受天线数目
    PRI = 4e-3;
    
    radar_range_profile = [];
    range_axis = linspace(0, fs * c / 2 / K, N);
    velo_axis = linspace(-lambda / 4 / Tc, lambda / 4 / Tc, M);
    time_axis = linspace(0, PRI * fnumber, n_chirps * n_samples);
    
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
    range_profile = [];
    for ttxx = 1:tx
        data_radar_1 = squeeze(t_lvds_data(rx_select_tx * rx - 3, :, :));   %RX1
        data_radar_2 = squeeze(t_lvds_data(rx_select_tx * rx - 2, :, :));   %RX2
        data_radar_3 = squeeze(t_lvds_data(rx_select_tx * rx - 1, :, :));   %RX3
        data_radar_4 = squeeze(t_lvds_data(rx_select_tx * rx, :, :));   %RX4    
            
        data_radar_1 = data_radar_1 - mean(data_radar_1, 2);
        data_radar_2 = data_radar_2 - mean(data_radar_2, 2);
        data_radar_3 = data_radar_3 - mean(data_radar_3, 2);
        data_radar_4 = data_radar_4 - mean(data_radar_4, 2);
        %% 3维FFT处理
        %距离FFT
        range_win = hamming(n_samples);   %加海明窗
        doppler_win = hamming(n_chirps);
        range_profile1 = abs(squeeze((fft(data_radar_1, N))));
        range_profile2 = abs(squeeze((fft(data_radar_2, N))));
        range_profile3 = abs(squeeze((fft(data_radar_3, N))));
        range_profile4 = abs(squeeze((fft(data_radar_4, N))));
        if isempty(range_profile)
            range_profile = range_profile1 + range_profile2 + range_profile3 + range_profile4;
        else
           range_profile = range_profile + range_profile1 + range_profile2 + range_profile3 + range_profile4;
        end
    end
    
    range_profile = range_profile / tx * rx;

        figure(10001)
        imagesc(abs(range_profile))
        %多普勒FFT
        radar_range_profile = [radar_range_profile, range_profile];
    end
 
    figure(10000)
    imagesc(time_axis, range_axis, db(radar_range_profile))
    xlabel("时间(s)")
    ylabel("距离(m)")
end