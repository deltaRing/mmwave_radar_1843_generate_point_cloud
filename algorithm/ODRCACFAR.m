% �������CACFAR
% ���룺1 map��RDmap RAmap ����
% ���룺2 R��    �ο���Ԫ
% ���룺3 P��    ������Ԫ
% ���룺4 gain�� ����
% �����1 detect_map��
% �����2 detect_result��
function [detect_map, detect_result] = ODRCACFAR(map, R, P)
    if nargin == 1
        R = 5;       % �ο���Ԫ����Ϊ 5
        P = 5;       % ������Ԫ����Ϊ 3
        Pfa = 0.2; % �龯�� 
    end
    
    size_x = size(map, 1);
    size_y = size(map, 2);
    L_slipper       = R + P; % �ܹ���Ҫ�ĵ�Ԫ
    L_slipper_P     = P;     % ������Ԫ
    detect_map      = zeros(size_x, size_y);    % �����ͼ
    detect_result   = [];    % �����
    
    for tt = 1:size_y
        for rr = 1:size_x
            test_unit_rr = rr + L_slipper;    % ��¼���Ե�Ԫ
            end_unit_rr = rr + 2 * L_slipper; % ��¼�����Ĳ��Ե�Ԫ
            if end_unit_rr > size_x, continue, end
            test_unit_ps = test_unit_rr - L_slipper_P;    % ��¼��ʼ��Ԫ�Լ�������Ԫ
            test_unit_pe = test_unit_rr + L_slipper_P;    % ��¼������Ԫ
            test_units = [map(rr:test_unit_ps, tt); ...
                map(test_unit_pe:end_unit_rr, tt)];     % ��¼��Ԫ
            if abs(map(test_unit_rr, tt)) > sum(abs(test_units)) * (Pfa^(-1/R/2)-1)
                detect_map(test_unit_rr, tt) = map(test_unit_rr, tt);
                detect_result = [detect_result; test_unit_rr, tt, map(test_unit_rr, tt)];
            end
        end
    end

end

