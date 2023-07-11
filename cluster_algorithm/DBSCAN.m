% DBSCAN 算法
% 输入1：数据集     data (N x 3)
% 输入2：聚类半径   EPS
% 输入3：最小数据量 minPoints
% 输出1：聚类数据   data_cluster
function data_cluster = DBSCAN(data, eps, minpoints)
    % 提供默认值
    if nargin == 1
        eps = 0.5;      % 最大半径是0.6
        minpoints = 25; % 最小点数是25个点
    end
    % 数据是否已经被聚类
    data_type    = zeros(size(data, 1), 1);
    data_points  = [];      % 聚类的记录点数
    % 数据标号
    data_index   = 1;
    for ii = 1:size(data, 1)
        % 初始化
        if data_type(ii) == 0
            data_type(ii)    = data_index;     % 数据类别
            points           = 1;              % 聚类的点数
            data_index       = data_index + 1; % 数据标号自增
        else
            points           = data_points(data_type(ii)); % 取出已有点数
        end
        for iii = 1:size(data, 1) 
            % 被聚过类了 或者是 就是它自己
            if ii == iii, continue, end
            if data_type(iii), continue, end
            if norm(data(iii, :) - data(ii, :)) < eps
                data_type(iii)    = data_type(ii);  % 类别=1
                points            = points + 1;     % 点数加1
            end
        end
        data_points(data_type(ii)) = points;        % 更新点数 
    end
    
    % 找到对应的数据标号，并删除对应的数据标号
    outliers_index              = find(data_points < minpoints);
    outliers_array              = [];
    for ii = 1:length(outliers_index)
        outliers_array = [outliers_array find(data_type == outliers_index(ii))'];
    end
    % 对数据进行擦除
    data(outliers_array,:)        = [];
    data_type(outliers_array,:)   = [];
    data_points(outliers_index)   = [];
    
    % 重新对数据进行编号 长官请讲
    re_index = 1;
    for ii = 1:data_index - 1
        index = find(data_type == ii);
        if ~isempty(index)
            data_type(index) = re_index;
            re_index = re_index + 1;
        end
    end
        
    data_cluster = [data data_type]; % 被聚类的数据
    
    if 0
        figure(10023)
        hold on
        for ii = 1:re_index - 1
            data_scatter = data_cluster(find(data_cluster(:,4)==ii), 1:3);
            scatter3(data_scatter(:,1),data_scatter(:,2),data_scatter(:,3))
        end
        title("DBSCAN聚类结果 EPS=0.5, point=25")
        xlabel("距离 (meter)")
        ylabel("速度 (m/s)")
        zlabel("角度 (rad)")
        legend
    end
end

