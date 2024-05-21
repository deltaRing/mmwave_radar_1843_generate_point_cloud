% MVDR 自适应波束形成算法
% 输入：发射天线所对应的距离像
% 输入：空间谱点数
% 输入：载频
% 输出：距离-角度谱图
function range_angle_map = MVDR(range_profile_tx, SpaceNum, f0)
    if nargin == 1
        SpaceNum = 512;      % 遍历次数
        f0       = 77 * 1e9; % 载频
    end
    warning('off')
    range_profile_tx = squeeze(range_profile_tx(:,1,:));    % 只读取一部分数据
    theta_axis       = linspace(-pi / 2, pi / 2, SpaceNum); % 空间谱
    tx_num           = size(range_profile_tx, 2);           % 发射天线数量
    c                = 3e8;                                 % 光速
    lambda           = c / f0;                              % 波长
    k                = 2 * pi * f0 / c;                     %
    space            = lambda / 2;                          % 天线间距
    P                = [1 : tx_num];                        %
    range_angle_map  = [];                                  % 返回
    for rr = 1:size(range_profile_tx, 1)
        yy = range_profile_tx(rr, :);
        R = yy' * yy; % 自相关矩阵
        for aa = 1:SpaceNum
            p = exp(1j*k*space*P*sin(theta_axis(aa)))';
%             Wcc = inv(R) * p / (p' / R * p);
%             B(aa) = Wcc' * R * Wcc;
            B(aa) = 1 / (p' * inv(R) * p);
        end
        range_angle_map = [range_angle_map; B];  % 计算功率谱
    end
    range_angle_map = db(range_angle_map);
end

