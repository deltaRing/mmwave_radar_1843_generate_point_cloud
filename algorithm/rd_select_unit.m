% ʹ��RD��ͼ��RA��ͼ����һ�����Ŀ������
% ���巽��������ʾ��
% ���ȶ�ÿ�����뵥Ԫ����CFAR ������ѡȡ���뵥Ԫ �ٶȵ�Ԫ �Լ� ��Ӧ�ĽǶȵ�Ԫ
% �������ϵ���ȡ��   
%   ��ÿ�����ݽ��У��޼ල������ ��RDA���Զ�����Ϊ
%   ����������� Ĭ�ϵ����ǵ���
%   ����Ϊ�������������˻����������ʼ��С�ڵ�ǽ��ľ��룩
% ����1��RD_MAP ����-��������ͼ Ӧ����һ������ N x RangeUnits x DopplerUnits
% ����2��RA_MAP ����-�Ƕ���ͼ Ӧ����һ������ N x RangeUnits x AngleUnits
% ����3��range_axis ������
% ����4��doppler_axis �ٶ���
% ����5��angle_axis �Ƕ���
% ����6��disable_velo ���ε��ٶ���
% �����
% ����Ȥ�ľ��뵥Ԫ��range_unit 
% ����Ȥ�ĽǶȵ�Ԫ��angle_unit 
function detect_units = rd_select_unit(rd_map,ra_map, ...
            range_axis,doppler_axis,angle_axis, ...
            drone_height, disable_velo)
    % 1.6��һ�µ��ٶ�
    if nargin == 5, disable_velo = 20; drone_height = 1.3; end 
    % ���˻��߶�����Ӧ�ľ��뵥Ԫ���
    drone_height_unit = sum(drone_height > range_axis) + 5; % ��ֹ�߶����������Լ������Ӱ��
    % �ҵ������ٶ�
    center_velo = length(doppler_axis) / 2; % �����ٶ�
    % ����Ȥ���ٶȵ�Ԫ
    valid_velo_min  = center_velo - disable_velo;
    valid_velo_max  = center_velo + disable_velo;
    % ѡȡ����Ȥ�ľ��뵥Ԫ���ٶȵ�Ԫ
    [detect_map_rd, detect_result_rd] = ODRCACFAR(rd_map);
    detect_threshold                  = mean(mean(detect_map_rd(find(detect_map_rd > 0))));
    detect_map_rd(find(detect_map_rd < detect_threshold)) = 0;
    % ������ص��ٶȵ�Ԫ
    detect_map_rd(:, 1:valid_velo_min)   = 0;
    detect_map_rd(:, valid_velo_max:end) = 0;
    % ������صľ��뵥Ԫ
    detect_map_rd(1:drone_height_unit, :) = 0;
    % �˳���صļ����
    detect_result_rd(find(detect_result_rd(:, 2) < valid_velo_min),:)    = [];  % �ٶ��˳�
    detect_result_rd(find(detect_result_rd(:, 2) > valid_velo_max),:)    = [];  % �ٶ��˳�
    detect_result_rd(find(detect_result_rd(:,3) < detect_threshold), :)  = [];  % ��ֵ�˳�
    detect_result_rd(find(detect_result_rd(:, 1) < drone_height_unit),:) = [];  % �����˳�
    % �������ϵļ����������RA��ͼ
    detect_result_ra = ra_select_unit(ra_map, detect_result_rd);
    % ������� ������ӳ�䵽����ľ�����
    angles = angle_axis(detect_result_ra);
    ranges = range_axis(detect_result_rd(:, 1));
    velos  = doppler_axis(detect_result_rd(:, 2));
    % �������
    detect_units = [ranges; velos; angles]';
end

