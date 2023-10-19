% ra选择距离单元
% 根据：遍历距离单元
%       选取最大的响应单元
% 输入1：RA_MAP        距离-角度谱图
% 输入2：Detect_result 检测结果 
% 输入3：Range_spread  遍历距离结果
% 输出： Angle_Result  检测结果
function angle_result = ra_select_unit(ra_map, detect_result, range_spread)
    if nargin == 2, range_spread = 5; end
    angle_result = [];
    for dd = 1:size(detect_result, 1)
        % 选择检测距离单元
        rr = detect_result(dd, 1);
        rr_min = rr - range_spread;
        rr_max = rr + range_spread;
        % 检测是否超出范围
        if rr_min < 1, rr_min = 1; end
        if rr_max > size(ra_map, 1), rr_max = size(ra_map, 1); end
        % 具体的距离单元
        ra_map_select = ra_map(rr_min:rr_max, :);
        [index_x, index_y] = find(max(max(abs(ra_map_select))) == abs(ra_map_select));
        angle_result = [angle_result index_y];
    end
end

