% 检查是否存在新的航迹
% 输入1：new_trace 新的航迹
% 输入2：previous_trace 旧的航迹结果
% 输出1：res 检测结果
% 输出2：新航迹的index 若有：不为空 若无：[] 空数组
% 输出3：被删除的index 若有：不为空 若无：[] 空数组
% 输出4：航迹丢失的index 若有：不为空 若无：[] 空数组 (只会在) 
function [res, new_index, remove_index, loss_index] = check_new_trace(new_trace, previous_trace)
    if isempty(new_trace) && isempty(previous_trace)
        res = 0;
        new_index = [];
        remove_index = [];
        loss_index = [];
        return
    end

    if ~isempty(new_trace) && isempty(previous_trace)
        res = 1;
        new_index = 1:size(new_trace, 1);
        remove_index = [];
        loss_index = [];
        return
    end

    if isempty(new_trace) && ~isempty(previous_trace)
        res = 1;
        new_index = [];
        remove_index = 1:size(previous_trace, 1);
        loss_index = [];
        return
    end

    res = 0; new_index = []; remove_index = []; loss_index = [];
    if size(previous_trace, 1) ~= size(new_trace, 1)
        res = 1;
    else
        
    end

    for ii = 1:size(previous_trace, 1)
        trace_id = previous_trace(ii, 5);
        % 先检查有没有被删的index
        res_exist = sum(find(new_trace(:, 5) == trace_id));
        if ~res_exist
            remove_index = [remove_index ii];
            res = 1;
        end
    end
        
    for ii = 1:size(new_trace, 1)
        trace_id = new_trace(ii, 5);
        % 再检查有没有新加入的index
        res_exist = sum(find(previous_trace(:, 5) == trace_id));
        if ~res_exist
            new_index = [new_index ii];
            res = 1;
        end
    end

    for ii = 1:size(new_trace, 1)
        trace_id = new_trace(ii, 5);
        % 最后检查有没有丢失航迹的index
        res_exist = find(previous_trace(:, 5) == trace_id);
        if ~isempty(res_exist)
            index_diff = new_trace(ii, 4) - previous_trace(res_exist, 4);
            if index_diff > 2
                loss_index = [loss_index ii];
            end
        end
    end
end