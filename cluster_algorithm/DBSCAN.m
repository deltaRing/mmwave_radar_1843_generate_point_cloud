% DBSCAN �㷨
% ����1�����ݼ�     data (N x 3)
% ����2������뾶   EPS
% ����3����С������ minPoints
% ���1����������   data_cluster
function data_cluster = DBSCAN(data, eps, minpoints)
    % �ṩĬ��ֵ
    if nargin == 1
        eps = 0.5;      % ���뾶��0.6
        minpoints = 25; % ��С������25����
    end
    % �����Ƿ��Ѿ�������
    data_type    = zeros(size(data, 1), 1);
    data_points  = [];      % ����ļ�¼����
    % ���ݱ��
    data_index   = 1;
    for ii = 1:size(data, 1)
        % ��ʼ��
        if data_type(ii) == 0
            data_type(ii)    = data_index;     % �������
            points           = 1;              % ����ĵ���
            data_index       = data_index + 1; % ���ݱ������
        else
            points           = data_points(data_type(ii)); % ȡ�����е���
        end
        for iii = 1:size(data, 1) 
            % ���۹����� ������ �������Լ�
            if ii == iii, continue, end
            if data_type(iii), continue, end
            if norm(data(iii, :) - data(ii, :)) < eps
                data_type(iii)    = data_type(ii);  % ���=1
                points            = points + 1;     % ������1
            end
        end
        data_points(data_type(ii)) = points;        % ���µ��� 
    end
    
    % �ҵ���Ӧ�����ݱ�ţ���ɾ����Ӧ�����ݱ��
    outliers_index              = find(data_points < minpoints);
    outliers_array              = [];
    for ii = 1:length(outliers_index)
        outliers_array = [outliers_array find(data_type == outliers_index(ii))'];
    end
    % �����ݽ��в���
    data(outliers_array,:)        = [];
    data_type(outliers_array,:)   = [];
    data_points(outliers_index)   = [];
    
    % ���¶����ݽ��б�� �����뽲
    re_index = 1;
    for ii = 1:data_index - 1
        index = find(data_type == ii);
        if ~isempty(index)
            data_type(index) = re_index;
            re_index = re_index + 1;
        end
    end
        
    data_cluster = [data data_type]; % �����������
    
    if 0
        figure(10023)
        hold on
        for ii = 1:re_index - 1
            data_scatter = data_cluster(find(data_cluster(:,4)==ii), 1:3);
            scatter3(data_scatter(:,1),data_scatter(:,2),data_scatter(:,3))
        end
        title("DBSCAN������ EPS=0.5, point=25")
        xlabel("���� (meter)")
        ylabel("�ٶ� (m/s)")
        zlabel("�Ƕ� (rad)")
        legend
    end
end

