% 使用XMEANS 进行聚类
% 输入1：等待聚类的数据 data
% 输入2：最小的聚类数目 Kmin
% 输入3：最大的聚类数目 Kmax
% 输出1：聚类中心 center
% 输出2：新的数据 data_new
function [center, data_new] = xmeans(data, kmax)
    % 单个数据的大小（n x 3 + 1）
    k = 1; % 从单个聚类开始迭代
    while 1
        ok = k;
        % 根据最小的Kmin来计算
        [center_init, data_init] = kmeans(data, k);
        index_num = size(data_init, 2);  % index所在的位置
        free_para = index_num - 1;       % 自由参数量
        
        o_bic_scores = [];
        for ii = 1:k
            % 分别计算 K = 1 以及 K = 2 的时候 BIC分数是多少
            data_select = data_init(find(data_init(:, index_num) == ii), 1:3);  % 选择的数据
            data_size   = size(data_init, 1);                                   % 数据总量
            center_ii   = center_init(ii,:);                                    % 到中心位置的点
            distance    = norm(data_select - center_ii)^2;                      % 到中心位置的距离
            score_1     = BIC_points(1, data_size, free_para, ...
                [distance], [data_size]);                                       % 计算BIC分数
            o_bic_scores  = [o_bic_scores score_1];
        end
        
        SK    = 2; % 计算子聚类
        ADD_K = 0; % 需要加入的K
        s_bic_scores = [];
        for ii = 1:k 
            data_select         = data_init(find(data_init(:, index_num) == ii), 1:3); % 选择的数据
            [center_2, data_2]  = kmeans(data_select, 2);
            data_select_1       = data_2(find(data_2(:, index_num) == 1), 1:3);   % 选择的数据 
            data_select_2       = data_2(find(data_2(:, index_num) == 2), 1:3);  
            data_size1          = size(data_select_1, 1);                              % 数据总量
            data_size2          = size(data_select_2, 1);                     
            distance1           = norm(data_select_1 - center_2(1, :))^2;              % 数据到各自中心的方差
            distance2           = norm(data_select_2 - center_2(2, :))^2;   
            score_2             = BIC_points(2, data_size, free_para, ...
                [distance1 distance2], [data_size1 data_size2]);    % 计算BIC分数
            s_bic_scores  = [s_bic_scores score_2];                 % 分数迭代部分
        end
        
        ADD_K = length(find(s_bic_scores > o_bic_scores));
        k = k + ADD_K;
        
        if ok == k || k >= kmax, break; end
    end
    % 计算聚类结果
    [center, data_new] = kmeans(data, k);
    
    if 0
        figure(10021)
        scatter3(data_select(:, 1),data_select(:, 2),data_select(:, 3))
        figure(10022)
        scatter3(data_select_1(:, 1),data_select_1(:, 2),data_select_1(:, 3))
        hold on
        scatter3(data_select_2(:, 1),data_select_2(:, 2),data_select_2(:, 3))
    end
end

% 计算BIC分数
% 输入1：K 簇的总数 （K = 1、2）
% 输入2：R 样本总数 （所有的点）
% 输入3：M 数据维度 （自由参数）(3)
% 输入4：distance 数据中心到各个数据的距离平方
% 输入5：cluster_size 每个集群的数量
% 输出1：score BIC分数 
function score = BIC_points(K, R, M, distance, cluster_size)
    var = 1 / (R - K) * distance;  % 方差
    L = 0;
    for ii = 1:K
        L = L + logLikelihood(K, R, cluster_size(ii), M, var(ii));
    end
    numParameters = (M + 1) * K; % BIC公式
    score = L - 0.5 * numParameters * log(R);
end

% 计算对数似然函数
% 输入1：K   簇的总数
% 输入2：R   样本总数
% 输入3：Ri  数据量
% 输入4：M   数据维度（3）
% 输入5：var 方差
% 输出1：logLike:对数函数
function loglike = logLikelihood(K, R, Ri, M, var)
    p1 = -Ri * log(2 * pi);
    p2 = -Ri * M * log(var);
    p3 = -(Ri - K);
    p4 = Ri * log(Ri);
    p5 = -Ri * log(R);
    loglike = (p1 + p2 + p3) / 2.0 + p4 + p5;
end
