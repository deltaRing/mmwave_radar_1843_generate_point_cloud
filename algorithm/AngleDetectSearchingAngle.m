%% cfar_detect_info CFAR �����Ϣ
% ���룺
%   1. �Ƕ�-������ͼ 
%   2. cfar�����
% �����
%   1. ��Ӧ�ĽǶ�
function detect_res = AngleDetectSearchingAngle(angle_map, ...
    cfar_detect_info, ...
    offset, ...
    disable_range)
if nargin < 3
   offset = 10; 
   disable_range = 60; % ���ǿɶ�̬������
end

% �����ǰ��disable��Ԫ��Ч����ô��
if disable_range > 0 && disable_range <= size(angle_map, 1)
    angle_map(1:disable_range, :) = 0;
end

detect_res = [];
% ���뵥ԪУ׼
cfar_detect_info = cfar_detect_info(:, 1) + offset;
% ������Щ�Ƕ�
angle_disable_min = 64;
angle_disable_max = size(angle_map, 2) - 64;
% ȫ�ֵķ���
global_mean_amp = mean(mean(angle_map(disable_range:end, :))); % ȫ�ֵ�ƽ������
global_max_amp = max(max(angle_map(disable_range:end, :)));    % ȫ�ֵ�������
if global_max_amp / global_mean_amp < 1.5     % ȷ����ǰ����ֵ֮���Ƿ�Ƚ�����ȷ̽�⵽Ŀ�� 
    return
end

    for ii = size(cfar_detect_info(:, 1)):-1:1
        % ɾ������̫�̵ĵ�
        if abs(cfar_detect_info(ii, 1) < disable_range)
            cfar_detect_info(ii, :) = [];
        end
    end

    cfar_detect_info = unique(cfar_detect_info(:, 1) + offset);
    for ii = 1:size(cfar_detect_info, 1)
        r_index = cfar_detect_info(ii);
        
        % Ѱ����Ӧ���ľ���
        r_index = range_search(angle_map, r_index, 3);
        
        angle_det_res = angle_map(r_index, :);
        %  get ������Ӧ
        max_amp = max(angle_det_res);
        %  get ƽ������Ӧ
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

%% Ѱ�����е�range�����ҵ����ķ��ȵľ��뵥Ԫ
function range_searched_index = range_search(az_map, range_index, offset)
    ranges_min = range_index - offset; ranges_max = range_index + offset;
    if ranges_min < 0 
        ranges_min = offset;
    end
    if ranges_max > size(az_map, 1)
        ranges_max = size(az_map, 1);
    end
    % �����趨
    ranges = [ranges_min:1:ranges_max];
    az_selected_map = az_map(ranges, :); % ѡ��ǶȾ�����
    detect_index = [];
    
    for ii = 1:length(ranges)
        [val, index] = max(az_map(ii, :)); % �ҵ����ֵ
        detect_index = [detect_index index];
    end
    
    % �����������Ϊ��ȡ�ú���ĽǶȵ�Ԫ���
    if length(unique(detect_index)) == length(detect_index)
        index = fix(mean(detect_index)); % ʹ��ƽ����
    else
        index = mode(detect_index); % ȡ������
    end
    
    range_searched_index = ranges(index);
end