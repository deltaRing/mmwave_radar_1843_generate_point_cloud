function list = DetectPointFilter(data)
    list = [];
    if isempty(data)
        fprintf("Data 需要三维的数据");
        return;
    end
    
    if length(size(data{1})) ~= 2
        fprintf("Data 需要三维的数据");
        return;
    end 
    
    max_size_data = 0;
    for rx = 1:length(data)
        if max_size_data < size(data{rx}, 1)
            max_size_data = size(data{rx}, 1);
        end
    end
    occ_matrix = zeros(length(data), max_size_data); % 占用矩阵 如果已经检测完毕了 那就置为1
    roi_list = []; % 感兴趣的距离单元
    voi_list = []; % 感兴趣的速度单元
    threshold_rr = 5; % 至少持续三个速度单元
    thresgate_rr = 3; % 在速度单元五格子以内都算是
    
    for rx = 1:length(data)
       % 每个天线
       for rr = 1:size(data{rx},1) % 每个检测结果
           if occ_matrix(rx, rr) > 0
               continue; % 检测过了
           end
           rrr = data{rx}(rr, 1);
           vvv = data{rx}(rr, 2);
           ddd = 1; % 检测次数
           occ_matrix(rx, rr) = 1; % 标记为已经遍历过了
           for rrxx = 1:length(data)
               for rrrr = 1:size(data{rrxx},1)
                   % 先看看是不是检测过了
                    if occ_matrix(rrxx, rrrr) > 0
                       continue; % 检测过了
                    end
                    % 如果速度小于速度门阈值
                    if abs(vvv - data{rrxx}(rrrr, 2)) <= thresgate_rr
                        occ_matrix(rrxx, rrrr) = 1; % 标记为1
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