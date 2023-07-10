% K-Means 聚类算法
% 输入：n x 3 维度的数据
% 输入：K 初始选择的数据
% 输出：聚类中心 Center
% 输出：聚类数据 Data n x 3 + 1 (x y z index)
function [center, data_new] = kmeans(data, K)
    % 随机选择初始化的数据
    init_index  = randi(size(data, 1), [1, K]);  % 随机的索引 
    init_data   = data(init_index, :);                       % 挑选的数据
    data_index  = [];                                        % 存储的数据加上索引
    for ii = 1:length(init_index), data_index = [data_index; init_data(ii, :) ii]; end
    data_size   = size(data_index, 2);                       % 数据大小
    
    while 1
        for ii = size(data, 1):-1:1
            data_distance = []; % 数据距离
            for jj = 1:K
                distance = norm(init_data(jj, :) - data(ii, :));  % 到K的距离
                data_distance = [data_distance distance];         % 存储数据
            end
            [value, index] = min(data_distance);            % 找到最近的距离
            data_index = [data_index; data(ii, :) index];   % 添加数据
        end
        
        % 重新计算聚类中心
        new_center = [];  % 新的聚类中心
        for iii = 1:K
            DOI = find(data_index(:, data_size) == iii);       % 寻找数据
            cluster_data = data_index(DOI, 1:data_size - 1);   % 选择数据
            new_center   = [new_center; mean(cluster_data)];
        end
        
        % 已经收敛了
        if norm(new_center - init_data) < 1e-6, break, end
        init_data = new_center; % 更新数据
        data_index  = [];       % 存储的数据加上索引
    end
    
    data_new = data_index;   % 存储的数据加上索引
    center   = new_center;   % 记录的数据中心
    
    % 绘制数据
    if 0
        figure(10020)
        hold on
        for ii = 1:K
            data_scatter = data_index(find(data_index(:,4)==ii), 1:3);
            scatter3(data_scatter(:,1),data_scatter(:,2),data_scatter(:,3))
        end
        title("K-means聚类结果（指定K=2）")
        xlabel("距离 (meter)")
        ylabel("速度 (m/s)")
        zlabel("角度 (rad)")
        legend
    end
end
