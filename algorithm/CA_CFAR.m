% 输入：
% map：距离-多普勒谱图
% R: 参考单元
% P：保护单元
% Pfa：虚警率
% 输出：
% detect_map：检测后的谱图
% detect_result：检测结果
function [detect_map, detect_result] = CA_CFAR(map, R, P, Pfa)
    if nargin == 1
        R = 8;       % 参考单元设置为 5
        P = 6;       % 保护单元设置为 3
        Pfa = 2e-2; % 虚警率 
    end

    size_x = size(map, 1);
    size_y = size(map, 2);
    L_slipper       = R + P; % 总共需要的单元
    L_slipper_P     = P;     % 保护单元
    detect_map      = zeros(size_x, size_y);    % 检测谱图
    detect_result   = [];    % 检测结果
    detect_threshold = [];
    compute_result = [];
    
    for tt = 1:size_y
        for rr = 1:size_x
            test_unit_rr = rr + L_slipper;    % 记录测试单元
            end_unit_rr = rr + 2 * L_slipper; % 记录结束的测试单元
            if end_unit_rr > size_x, continue, end
            test_unit_ps = test_unit_rr - L_slipper_P;    % 记录起始单元以及结束单元
            test_unit_pe = test_unit_rr + L_slipper_P;    % 记录结束单元
            test_units = [map(rr:test_unit_ps, tt); ...
                map(test_unit_pe:end_unit_rr, tt)];     % 记录单元
            if abs(map(test_unit_rr, tt)) > sum(abs(test_units)) * (Pfa^(-1/R/2)-1)
                detect_map(test_unit_rr, tt) = map(test_unit_rr, tt);
                detect_result = [detect_result; test_unit_rr, tt, map(test_unit_rr, tt)];
            end
            detect_threshold(test_unit_rr, tt) = sum(abs(test_units)) * (Pfa^(-1/R/2)-1);
            compute_result(test_unit_rr, tt) = abs(map(test_unit_rr, tt));
        end
    end

end