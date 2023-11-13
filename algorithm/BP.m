% BP 成像
% 论文来源：https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553059
% 输入1：所有的距离像 antennas x ADC_Samples
% 输入2：距离尺度  range_axis 1 x FFT
% 输入3：波长 lambda 1x1
% 输入4：天线位置 antennas txs * rxs x 2
% 输入5：txs 发射天线数目 1 x 1
% 输入6：rxs 接收天线数目 1 x 1
% 输入7：img_size 成像图大小 1 x 1
% 输入8：range_max 最大距离 1 x 1 
function img = BP(range_profiles, range_axis, lambda, antennas, txs, rxs, img_size, range_max)
    img = zeros(img_size, img_size); % 对应距离？ 0---10m
    range_axis_x = linspace(-range_max / 2, range_max / 2, img_size);
    range_axis_y = linspace(0, range_max, img_size);
    channels = zeros(img_size, img_size, txs * rxs);
    for ii = 1:txs
        for jj = 1:rxs
            ant = antennas((ii - 1) * rxs + jj, :); % 天线位置

            % ------------------> 方位向
            % | []
            % |   |
            % |  \|/   <-- 计算当前点tau（延迟时间） 叠加range_profile
            % | [] ← 遍历下一个单元 （最后遍历下一个方位向）  
            % |
            % |
            % |
            % |
            % \/ 距离向

            for iii = 1:img_size
                for jjj = 1:img_size
                     range_cell = norm([range_axis_x(jjj) range_axis_y(iii)] - ant);
                     kkk = abs(range_cell - range_axis); % 找到最接近值
                     kkk_min = min(kkk);
                     kkk = find(kkk == kkk_min);
                     if length(kkk) > 1, kkk = kkk(1); end
                     range_profile = range_profiles(kkk, (ii - 1) * rxs + jj);
                     img(iii,jjj) = img(iii, jjj) + range_profile * exp(1j * 4 * pi * range_cell / lambda);
                     channels(iii, jjj, ii * rxs + jj) = range_profile * exp(1j * 4 * pi * range_cell / lambda);
                end
            end
        end
    end

    % pcf
    data_PCF_sign = sign(channels);
    sigma_fai = std(data_PCF_sign,1,3);
    EPCF = (1-sigma_fai).^2;
    img = img.*EPCF;
end
