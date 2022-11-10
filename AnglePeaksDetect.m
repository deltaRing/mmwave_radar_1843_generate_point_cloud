function angle_detect = AnglePeaksDetect(angle_range_map)
    angle_detect = [];
    if length(size(angle_range_map)) ~= 2
        fprintf("需要 距离-角度2维谱");
        return;
    end
    
    max_peak = max(max(abs(angle_range_map)));
    mean_peak = mean(mean(abs(angle_range_map)));
    if mean_peak > max_peak * 0.5
        fprintf("警告 角度谱可能存在问题");
        return;
    end
    threshold_peak = max_peak * 0.3;
    for rr = 1:size(angle_range_map, 1)
        [peaks, aa] = findpeaks(abs(angle_range_map(rr, :)));
        for pp = 1:length(peaks)
            if peaks(pp) == 1 || peaks(pp) == size(angle_range_map, 2)
                continue;
            end
            if peaks(pp) >= threshold_peak
                angle_detect = [angle_detect; rr, aa(pp)];
            end
        end
    end
end