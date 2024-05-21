% 利用距离-角度谱得到角度信息
% 输入1：凝聚后检测结果
% 输入2：距离-角度谱图
% 输出1：距离信息
% 输出2：角度信息
function [range, angle] = getAngleInfo(res, azmap)
    range = []; angle = {};
    angle_FFT = size(azmap, 2);
    range_FFT = size(azmap, 1);
    
    if isempty(res), return, end
    res_range = fix(res(:, 1)); % 得到距离信息

    for ii = 1:length(res_range)
        azmap_data = azmap(ii, :).';
        [angle_value, angle_index] = findpeaks(azmap_data);
        % figure(1)
        % plot(azmap_data)
        % 角度处理
        % mean_value = mean(angle_value);
        % remove_index = find(mean_value > angle_value);
        % angle_index(remove_index) = [];
        % angle_value(remove_index) = [];
        angle_max_index = find(angle_value == max(angle_value));
        angle_index = angle_index(angle_max_index);
        range = [range res_range(ii)];
        if isempty(angle), angle = {angle_index};
        else, angle{length(angle) + 1} = angle_index; end
    end
end