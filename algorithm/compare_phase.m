% 利用比相法计算目标角度
% 输入1： rd_maps  range_index x velo_index x 通道数
% 输入2： cfar_res cacfar后的结果 range_index x velo_index
% 可选输入1： mti_offset mti的偏移量 
% 可选输入2： lambda 波长
% 可选输入3： d 天线间距 
% 输出1： 确定目标的距离
% 输出2： 确定目标的角度
function [range, angle] = compare_phase(rd_maps, cfar_res, mti_offset, lambda, d)
if nargin == 2
    mti_offset=10;
    lambda=0.0039;      %雷达信号波长
    d=lambda/2;       %天线阵列间距
end

% 初始化距离、角度
range = [];
angle = [];

% 距离单元校准
cfar_res(:, 1) = cfar_res(:, 1) + mti_offset;
% CFAR 遍历
for iii = 1:size(cfar_res, 1)
% 通道相位
    phase1 = phase(rd_maps(cfar_res(iii, 1), cfar_res(iii, 2), 1));
    phase2 = phase(rd_maps(cfar_res(iii, 1), cfar_res(iii, 2), 2));
    phase3 = phase(rd_maps(cfar_res(iii, 1), cfar_res(iii, 2), 3));
    phase4 = phase(rd_maps(cfar_res(iii, 1), cfar_res(iii, 2), 4));

    if phase1 < 0 || phase2 < 0 || phase3 < 0 || phase4 < 0
        continue;
    end

    thi1 = asin((phase2 - phase1)* lambda /(2*pi*1*d));
    thi2 = asin((phase3 - phase2)* lambda /(2*pi*1*d));
    thi3 = asin((phase4 - phase3)* lambda /(2*pi*1*d));
    thi4 = asin((phase4 - phase2)* lambda /(2*pi*2*d));

    angle_mean = (thi1 + thi2 + thi3 + thi4) / 4;
    angle = [angle angle_mean];
    range = [range cfar_res(iii, 1)];
end

end
