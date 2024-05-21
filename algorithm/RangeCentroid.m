% 距离凝聚函数
% 输入1：DetectRes（rr tt amp）
% 输入2：minimal_confirm 最小确认距离点数
% 输入3：consequence_range 连续距离半径
% 输入4：consequence_velo 连续速度半径
% 输入5：是否
% 输出1：Res 得到的距离凝聚结果 （rr tt var_rr var_tt amp）
%
function Res = RangeCentroid(DetectRes, ...
    minimal_confirm, ...
    consequence_range, ...
    consequence_velo)
    if nargin == 1
        minimal_confirm = 3;
        consequence_range = 3;
        consequence_velo = 3;
    end
    Res = [];

    while ~isempty(DetectRes)
        compare_range  = DetectRes(1, 1);      % 当前质心
        compare_velo   = DetectRes(1, 2);      % 当前速度
        compare_weight = abs(DetectRes(1, 3)); % 当前幅度
        temporary_compare_range  = DetectRes(1, 1);      % 当前质心
        temporary_compare_velo   = DetectRes(1, 2);      % 当前速度
        temporary_compare_weight = abs(DetectRes(1, 3)); % 当前幅度
        DetectRes(1, :) = [];            % 删除该数据
        delete_index = [];               % 有效的数据
        record_num = 1;                  % 记录观测次数
        Temporary_DetectRes = DetectRes; % 临时记录表
        % 融合周边目标 只有满足 当前质心
        while 1
            find_new_target = 0;
            temporary_delete_index = [];
            weight_record = [temporary_compare_weight];
            index = 1;
            while index < size(Temporary_DetectRes, 1) && size(Temporary_DetectRes, 1) >= 2
                cur_range  = Temporary_DetectRes(index, 1);      % 当前质心
                cur_velo   = Temporary_DetectRes(index, 2);      % 当前速度
                cur_weight = abs(Temporary_DetectRes(index, 3)); % 当前幅度 

                if abs(cur_range - temporary_compare_range) < consequence_range && ...
                    abs(cur_velo - temporary_compare_velo) < consequence_velo
                    find_new_target = 1;
                    temporary_delete_index = [temporary_delete_index index];
                    weight_record = [weight_record cur_weight];
                    record_num = record_num + 1;
                end
                index = index + 1;
            end

            if ~isempty(temporary_delete_index)
                temporary_compare_weight = sum(weight_record) / length(weight_record);
                weight_record = weight_record / sum(weight_record);
                temporary_compare_range = sum([temporary_compare_range Temporary_DetectRes(temporary_delete_index, 1)'] .* weight_record);
                temporary_compare_velo  = sum([temporary_compare_velo Temporary_DetectRes(temporary_delete_index, 2)'] .* weight_record);
                Temporary_DetectRes(temporary_delete_index, :) = [];
                delete_index = [delete_index temporary_delete_index];
                temporary_delete_index = [];
                weight_record = [temporary_compare_weight];
            end

            if ~find_new_target
                break;
            end
        end
        
        if record_num >= minimal_confirm
            range  = [compare_range, abs(DetectRes(delete_index, 1))'];
            velo   = [compare_velo, abs(DetectRes(delete_index, 2))'];
            weight = [compare_weight, abs(DetectRes(delete_index, 3))'];
            compare_weight = sum(weight) / length(compare_weight);
            weight = weight / sum(weight);
            compare_range = sum(range .* weight);
            compare_velo = sum(velo .* weight);
            var_range = var(range- compare_range);
            var_velo = var(velo- compare_velo);
            DetectRes(delete_index, :) = [];
            Res = [Res; compare_range compare_velo var_range, var_velo, compare_weight];
        end
    end
end