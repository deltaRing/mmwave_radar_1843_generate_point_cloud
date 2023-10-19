% raѡ����뵥Ԫ
% ���ݣ��������뵥Ԫ
%       ѡȡ������Ӧ��Ԫ
% ����1��RA_MAP        ����-�Ƕ���ͼ
% ����2��Detect_result ����� 
% ����3��Range_spread  ����������
% ����� Angle_Result  �����
function angle_result = ra_select_unit(ra_map, detect_result, range_spread)
    if nargin == 2, range_spread = 5; end
    angle_result = [];
    for dd = 1:size(detect_result, 1)
        % ѡ������뵥Ԫ
        rr = detect_result(dd, 1);
        rr_min = rr - range_spread;
        rr_max = rr + range_spread;
        % ����Ƿ񳬳���Χ
        if rr_min < 1, rr_min = 1; end
        if rr_max > size(ra_map, 1), rr_max = size(ra_map, 1); end
        % ����ľ��뵥Ԫ
        ra_map_select = ra_map(rr_min:rr_max, :);
        [index_x, index_y] = find(max(max(abs(ra_map_select))) == abs(ra_map_select));
        angle_result = [angle_result index_y];
    end
end

