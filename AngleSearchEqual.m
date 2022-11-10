function merged_res = AngleSearchEqual(angle_det_az, angle_det_el)
    merged_res = [];
    if isempty(angle_det_az) || isempty(angle_det_el)
        fprintf("��Ҫ��ά������\n");
        return;
    end
    
    if length(size(angle_det_az)) ~= 2 || length(size(angle_det_el)) ~= 2
        fprintf("��Ҫ��ά������\n");
        return;
    end

    angle_has_merged = zeros(1, size(angle_det_el, 1));
    for ii = 1:size(angle_det_az, 1)
        rr = angle_det_az(ii, 1);
        for jj = 1:size(angle_det_el, 1)
            if rr == angle_det_el(jj, 1)
                if angle_has_merged(jj)
                   continue;
                else
                    merged_res = [merged_res; rr, angle_det_az(ii, 2),...
                        angle_det_el(jj, 2)];
                    angle_has_merged(jj) = 1;
                    break;
                end
            elseif rr < angle_det_el(jj, 1)
                break;
            end
        end
    end

end