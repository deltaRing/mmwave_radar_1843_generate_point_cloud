% 速度向的CACFAR
% 输入：1 map：RDmap RAmap 均可
% 输入：2 R：    参考单元
% 输入：3 P：    保护单元
% 输入：4 gain： 增益
% 输出：1 detect_map：
% 输出：2 detect_result：
function [detect_map, detect_result] = ODVCACFAR(map, R, P, Pfa)
    if nargin == 1
        R = 5;     % 参考单元设置为 5
        P = 5;     % 保护单元设置为 3
        Pfa = 0.01; % 虚警率 
    end
    
    size_x = size(map, 1);
    size_y = size(map, 2);
    L_slipper       = R + P; % 总共需要的单元
    L_slipper_P     = P;     % 保护单元
    detect_map      = zeros(size_x, size_y);    % 检测谱图
    detect_result   = [];    % 检测结果
    
    for rr = 1:size_x
        for tt = 1:size_y
            test_unit_tt = tt + L_slipper;    % 记录测试单元
            end_unit_tt = tt + 2 * L_slipper; % 记录结束的测试单元
            if end_unit_tt > size_y, continue, end
            test_unit_ps = test_unit_tt - L_slipper_P;    % 记录起始单元以及结束单元
            test_unit_pe = test_unit_tt + L_slipper_P;    % 记录结束单元
            test_units = [map(rr, tt:test_unit_ps) map(rr, test_unit_pe:end_unit_tt)];     % 记录单元
            if abs(map(rr, test_unit_tt)) > sum(abs(test_units)) * (Pfa^(-1/R/2)-1)
                detect_map(rr, test_unit_tt) = map(rr, test_unit_tt);
                detect_result = [detect_result; rr test_unit_tt, map(rr, test_unit_tt)];
            end
        end
    end

end

