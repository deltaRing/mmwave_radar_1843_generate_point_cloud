function list = DetectPointFilter(data)
    list = [];
    if isempty(data)
        fprintf("Data ��Ҫ��ά������");
        return;
    end
    
    if length(size(data{1})) ~= 2
        fprintf("Data ��Ҫ��ά������");
        return;
    end 
    
    max_size_data = 0;
    for rx = 1:length(data)
        if max_size_data < size(data{rx}, 1)
            max_size_data = size(data{rx}, 1);
        end
    end
    occ_matrix = zeros(length(data), max_size_data); % ռ�þ��� ����Ѿ��������� �Ǿ���Ϊ1
    roi_list = []; % ����Ȥ�ľ��뵥Ԫ
    voi_list = []; % ����Ȥ���ٶȵ�Ԫ
    threshold_rr = 5; % ���ٳ��������ٶȵ�Ԫ
    thresgate_rr = 3; % ���ٶȵ�Ԫ��������ڶ�����
    
    for rx = 1:length(data)
       % ÿ������
       for rr = 1:size(data{rx},1) % ÿ�������
           if occ_matrix(rx, rr) > 0
               continue; % ������
           end
           rrr = data{rx}(rr, 1);
           vvv = data{rx}(rr, 2);
           ddd = 1; % ������
           occ_matrix(rx, rr) = 1; % ���Ϊ�Ѿ���������
           for rrxx = 1:length(data)
               for rrrr = 1:size(data{rrxx},1)
                   % �ȿ����ǲ��Ǽ�����
                    if occ_matrix(rrxx, rrrr) > 0
                       continue; % ������
                    end
                    % ����ٶ�С���ٶ�����ֵ
                    if abs(vvv - data{rrxx}(rrrr, 2)) <= thresgate_rr
                        occ_matrix(rrxx, rrrr) = 1; % ���Ϊ1
                        ddd = ddd + 1;
                    end
               end
           end
           if ddd >= threshold_rr
                roi_list = [roi_list rrr];
                voi_list = [voi_list vvv];
           end
       end
    end
    
    list = [roi_list; voi_list];
end