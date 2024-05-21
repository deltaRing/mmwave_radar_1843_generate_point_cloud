% MVDR ����Ӧ�����γ��㷨
% ���룺������������Ӧ�ľ�����
% ���룺�ռ��׵���
% ���룺��Ƶ
% ���������-�Ƕ���ͼ
function range_angle_map = MVDR(range_profile_tx, SpaceNum, f0)
    if nargin == 1
        SpaceNum = 512;      % ��������
        f0       = 77 * 1e9; % ��Ƶ
    end
    warning('off')
    range_profile_tx = squeeze(range_profile_tx(:,1,:));    % ֻ��ȡһ��������
    theta_axis       = linspace(-pi / 2, pi / 2, SpaceNum); % �ռ���
    tx_num           = size(range_profile_tx, 2);           % ������������
    c                = 3e8;                                 % ����
    lambda           = c / f0;                              % ����
    k                = 2 * pi * f0 / c;                     %
    space            = lambda / 2;                          % ���߼��
    P                = [1 : tx_num];                        %
    range_angle_map  = [];                                  % ����
    for rr = 1:size(range_profile_tx, 1)
        yy = range_profile_tx(rr, :);
        R = yy' * yy; % ����ؾ���
        for aa = 1:SpaceNum
            p = exp(1j*k*space*P*sin(theta_axis(aa)))';
%             Wcc = inv(R) * p / (p' / R * p);
%             B(aa) = Wcc' * R * Wcc;
            B(aa) = 1 / (p' * inv(R) * p);
        end
        range_angle_map = [range_angle_map; B];  % ���㹦����
    end
    range_angle_map = db(range_angle_map);
end

