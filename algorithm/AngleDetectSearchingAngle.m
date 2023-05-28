%% cfar_detect_info CFAR 检测信息
% 输入：
%   1. 角度-距离谱图 
%   2. cfar检测结果
% 输出：
%   1. 对应的角度
function detect_res = AngleDetectSearchingAngle(angle_map, ...
    cfar_detect_info, ...
    offset, ...
    disable_range)
if nargin < 3
   offset = 10; 
   disable_range = 60; % 这是可动态调整的
end

% 如果当前的disable单元有效，那么：
if disable_range > 0 && disable_range <= size(angle_map, 1)
    angle_map(1:disable_range, :) = 0;
end

detect_res = [];
% 距离单元校准
cfar_detect_info = cfar_detect_info(:, 1) + offset;
% 舍弃那些角度
angle_disable_min = 64;
angle_disable_max = size(angle_map, 2) - 64;
% 全局的幅度
global_mean_amp = mean(mean(angle_map(disable_range:end, :))); % 全局的平均幅度
global_max_amp = max(max(angle_map(disable_range:end, :)));    % 全局的最大幅度
if global_max_amp / global_mean_amp < 1.5     % 确定当前幅度值之差是否比较能明确探测到目标 
    return
end

    for ii = size(cfar_detect_info(:, 1)):-1:1
        % 删除距离太短的点
        if abs(cfar_detect_info(ii, 1) < disable_range)
            cfar_detect_info(ii, :) = [];
        end
    end

    cfar_detect_info = unique(cfar_detect_info(:, 1) + offset);
    for ii = 1:size(cfar_detect_info, 1)
        r_index = cfar_detect_info(ii);
        
        % 寻找响应最大的距离
        r_index = range_search(angle_map, r_index, 3);
        
        angle_det_res = angle_map(r_index, :);
        %  get 最大的响应
        max_amp = max(angle_det_res);
        %  get 平均的响应
        mean_amp = mean(angle_det_res);
        [a_values, a_index] = findpeaks(angle_det_res);
        valid_index = find(a_values > max_amp * 0.9 ...
            & a_values > global_mean_amp ...
            & mean_amp > global_mean_amp * 0.8);
        a_index = a_index(valid_index);
        
        for aa = 1:size(a_index, 2)
            if a_index(aa) < angle_disable_min || a_index(aa) > angle_disable_max
                continue
            end
            detect_res = [detect_res; r_index a_index(aa)];
        end
    end
end

%% 寻找所有的range，并找到最大的幅度的距离单元
function range_searched_index = range_search(az_map, range_index, offset)
    ranges_min = range_index - offset; ranges_max = range_index + offset;
    if ranges_min < 0 
        ranges_min = offset;
    end
    if ranges_max > size(az_map, 1)
        ranges_max = size(az_map, 1);
    end
    % 距离设定
    ranges = [ranges_min:1:ranges_max];
    az_selected_map = az_map(ranges, :); % 选择角度距离谱
    detect_index = [];
    
    for ii = 1:length(ranges)
        [val, index] = max(az_map(ii, :)); % 找到最大值
        detect_index = [detect_index index];
    end
    
    % 这个函数就是为了取得合理的角度单元标号
    if length(unique(detect_index)) == length(detect_index)
        index = fix(mean(detect_index)); % 使用平均数
    else
        index = mode(detect_index); % 取得众数
    end
    
    range_searched_index = ranges(index);
end