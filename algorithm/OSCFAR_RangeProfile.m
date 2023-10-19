% 针对距离像的OS-CFAR算法
% 输入：距离像 RangeProfile
% 输入：cell_num 网格单元值
% 输入：gain 增益值
% 输入：L：最低可接受单元
% 输入：H：最高可接受单元
% 输出：result：检测结果
function result = OSCFAR_RangeProfile(rangeprofile, cell_num, gain, L, H)
    if nargin < 2
        cell_num = 32;
        gain = 0.075; % 随便设置的
        L = 8;
        H = 24;
    end
    
    result = []; % 结果显示
    for tt = 1:size(rangeprofile, 2)
        for rr = 1:size(rangeprofile, 1) 
            if rr + cell_num > size(rangeprofile, 1)
                break; % 检测数组是否超标，这是一个距离上的OS-CFAR
            end
            test_index       = rr + cell_num / 2;                 % 测试序号 
            test_units       = rangeprofile(test_index, tt);      % 测试单元
            back_ground_unit = [rangeprofile(rr:test_index - 1, tt) ...
                rangeprofile(test_index+1:rr + cell_num, tt)];    % 背景单元
            sorted_units     = sort(abs(back_ground_unit));       % 排序后单元
            selected_units   = sorted_units(L:H);                 % 选择后单元
            sum_units        = sum(selected_units);               % 求和
            gain_units       = gain * sum_units;                  % 叠加增益
            if abs(test_units) >= gain_units
                result(rr, tt) = test_units;
            else
                result(rr, tt) = 0;
            end
        end
    end
    
end