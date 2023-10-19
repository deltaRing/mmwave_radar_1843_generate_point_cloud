% 使用RD谱图、RA谱图联合一起估计目标区域
% 具体方法如下所示：
% 首先对每个距离单元进行CFAR 初步地选取距离单元 速度单元 以及 对应的角度单元
% 经过以上的提取：   
%   对每个数据进行（无监督）聚类 （RDA）自动划分为
%   对最近的数据 默认当作是地面
%   （因为这里的情况是无人机到地面距离始终小于到墙面的距离）
% 输入1：RD_MAP 距离-多普勒谱图 应该是一个向量 N x RangeUnits x DopplerUnits
% 输入2：RA_MAP 距离-角度谱图 应该是一个向量 N x RangeUnits x AngleUnits
% 输入3：range_axis 距离轴
% 输入4：doppler_axis 速度轴
% 输入5：angle_axis 角度轴
% 输入6：disable_velo 屏蔽的速度轴
% 输出：
% 感兴趣的距离单元：range_unit 
% 感兴趣的角度单元：angle_unit 
function detect_units = rd_select_unit(rd_map,ra_map, ...
            range_axis,doppler_axis,angle_axis, ...
            drone_height, disable_velo)
    % 1.6米一下的速度
    if nargin == 5, disable_velo = 20; drone_height = 1.3; end 
    % 无人机高度所对应的距离单元求解
    drone_height_unit = sum(drone_height > range_axis) + 5; % 防止高度噪声带来对检测结果的影响
    % 找到中心速度
    center_velo = length(doppler_axis) / 2; % 中心速度
    % 感兴趣的速度单元
    valid_velo_min  = center_velo - disable_velo;
    valid_velo_max  = center_velo + disable_velo;
    % 选取感兴趣的距离单元、速度单元
    [detect_map_rd, detect_result_rd] = ODRCACFAR(rd_map);
    detect_threshold                  = mean(mean(detect_map_rd(find(detect_map_rd > 0))));
    detect_map_rd(find(detect_map_rd < detect_threshold)) = 0;
    % 屏蔽相关的速度单元
    detect_map_rd(:, 1:valid_velo_min)   = 0;
    detect_map_rd(:, valid_velo_max:end) = 0;
    % 屏蔽相关的距离单元
    detect_map_rd(1:drone_height_unit, :) = 0;
    % 滤除相关的检测结果
    detect_result_rd(find(detect_result_rd(:, 2) < valid_velo_min),:)    = [];  % 速度滤除
    detect_result_rd(find(detect_result_rd(:, 2) > valid_velo_max),:)    = [];  % 速度滤除
    detect_result_rd(find(detect_result_rd(:,3) < detect_threshold), :)  = [];  % 幅值滤除
    detect_result_rd(find(detect_result_rd(:, 1) < drone_height_unit),:) = [];  % 距离滤除
    % 根据以上的检测结果，遍历RA谱图
    detect_result_ra = ra_select_unit(ra_map, detect_result_rd);
    % 组合数据 将其重映射到具体的距离轴
    angles = angle_axis(detect_result_ra);
    ranges = range_axis(detect_result_rd(:, 1));
    velos  = doppler_axis(detect_result_rd(:, 2));
    % 打包数据
    detect_units = [ranges; velos; angles]';
end

