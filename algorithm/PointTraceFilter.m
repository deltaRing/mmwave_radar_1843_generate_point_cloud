% 点迹滤波器 （本算法适用于滤除雷达多余的杂点以及无规律的噪声点）
% 输入：帧id 新探测结果 确认探测结果 临时探测结果（记忆性结果） 近邻滤波阈值 添加次数
% 输出：新的确认探测结果[N x 4] 临时探测结果[N x 5]
% 假定所有的探测结果是：N x 3 的数据结构 （新的探测结果 确认探测结果的组成）
%    所有的临时探测结果是：N x 5 的数据结构：xyz 观测数目 最后观测帧的id
function [detect_res, new_temp_detect_res] = PointTraceFilter(frame_id, ...
    new_detect_res, ... % 新的探测结果 N x 3
    confirm_detect_res, ... % 确认的探测结果 N x 4
    temp_detect_res, ... % 临时的探测结果 N x 5
    same_point_range, ... % 近邻滤波阈值半径 single number
    minimal_number, ...  % 最小确认数目 single number
    expired_number)  % 临时存储结果过期数 single number
if nargin <= 4
    same_point_range = 0.65; % 0.65米以内的算作是一个目标
    minimal_number = 5;    % 最起码需要三次确认
    expired_number = 20;   % 如果连续二十帧没有出现 去除 temp_detect_res 和 confirm_detect_res 的结果
end

% 先融合所有的数据 （针对新探测结果）
% 按照融合半径：same_point_range
delete_list = [];
confirm_list = zeros(1, size(new_detect_res, 1));
for ii = size(new_detect_res, 1):-1:1
    if ii == 1
        break % 如果数组只有一个目标 直接跳过
    end
    if confirm_list(ii)
        continue % 如果确认删除 跳过
    end
    current_tar = new_detect_res(ii, :); % 先获取该点
    for jj = ii-1:-1:1
        if confirm_list(jj)
            continue;
        end
        temp_tar = new_detect_res(ii, :); % 计算距离
        range_of_new_tar = abs(norm(current_tar - temp_tar)); % 计算距离
        if range_of_new_tar < same_point_range
            delete_list = [delete_list jj]; % 删除列表
            confirm_list(jj) = 1;
        end
    end
end

clear confirm_list;
if ~isempty(delete_list)
    new_detect_res(delete_list, :) = []; % 删除数据
end

delete_list = [];
% 滤波后：去除过期的临时探测结果
for ii = size(temp_detect_res, 1):-1:1
    if abs(temp_detect_res(ii, 5) - frame_id) >= expired_number
        delete_list = [delete_list ii];
    end
end
if ~isempty(delete_list)
    temp_detect_res(delete_list, :) = []; % 删除数据
end

% 对所有数据进行近邻滤波
% 遍历temp_detect_res
if size(temp_detect_res, 1) > 0
    for ii = size(temp_detect_res, 1):-1:1
        data_prior = temp_detect_res(ii, 1:3);
        data_index = []; % 用来记录滤波的结果
        for jj = size(new_detect_res, 1):-1:1
            data_iter = new_detect_res(jj, :);
            if abs(norm(data_iter - data_prior)) > same_point_range
                continue; % 这个点不是邻近点
            end
            data_index = [data_index jj];
        end
        if ~isempty(data_index)
            for iii = 1:length(data_index)
                temp_detect_res(ii, 1:3) = temp_detect_res(ii, 1:3) + ...
                    new_detect_res(data_index(iii), :);
            end
            temp_detect_res(ii, 1:3) = temp_detect_res(ii, 1:3) / ...
                (length(data_index) + 1);
            temp_detect_res(ii, 4) = temp_detect_res(ii, 4) + 1; % 数据 + 1
            temp_detect_res(ii, 5) = frame_id;                   % 更新帧 id
            new_detect_res(data_index, :) = []; % 舍弃这些数据
        end
    end   
else
    % 探测结果
    for ii = size(new_detect_res, 1):-1:1
        temp_detect_res = [temp_detect_res; new_detect_res(ii, :) 1 frame_id];
        new_detect_res(ii, :) = []; % 删除数据
    end
end

% 滤除confirm_detect_res的结果
delete_list = [];
for ii = 1:size(confirm_detect_res, 1)
    if abs(confirm_detect_res(ii, 4) - frame_id) >= expired_number
        delete_list = [delete_list ii];
    end
end
confirm_detect_res(delete_list, :) = [];
    
% 滤波完成 先看看剩余的new_detect是否还能与确认航迹关联上
for ii = 1:size(confirm_detect_res, 1)
    conf_detect_res = confirm_detect_res(ii, 1:3);
    data_index = [];
    for jj = size(new_detect_res, 1):-1:1
        curr_detect_res = new_detect_res(jj, :);
        if abs(norm(conf_detect_res - curr_detect_res)) > same_point_range
            continue; % 太远了 不是可能的航迹点
        end
        data_index = [data_index jj];
        confirm_detect_res(ii, 1:3) = new_detect_res(jj, :); % 更新目标位置
        confirm_detect_res(ii, 4) = frame_id;                % 帧ID
        break;
    end
    if ~isempty(data_index)
        new_detect_res(data_index, :) = [];
    end
end
    
% 检查temp_detect_res 确认是否能加入confirm_detect_res 中
delete_list = [];
for ii = 1:size(temp_detect_res, 1)
    if temp_detect_res(ii, 4) >= minimal_number
        delete_list = [delete_list ii];
        temp_data = [temp_detect_res(ii, 1:3) frame_id];
        confirm_detect_res = [confirm_detect_res; temp_data];
    end
end
temp_detect_res(delete_list, :) = []; % 删除数据

% 把剩下的new 检测数据加入到temp_detect列表里面
for ii = size(new_detect_res, 1):-1:1
    temp_detect_res = [temp_detect_res; new_detect_res(ii, :) 1 frame_id];
    new_detect_res(ii, :) = [];
end

detect_res = confirm_detect_res;
new_temp_detect_res = temp_detect_res;
    
end
