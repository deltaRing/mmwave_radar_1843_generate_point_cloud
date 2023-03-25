function angle_detect = AnglePeaksDetect(angle_range_map)
    angle_detect = [];
    range_index = [];
    if length(size(angle_range_map)) ~= 2
        fprintf("需要 距离-角度2维谱");
        return;
    end
    
    ratio = 0.85;
    max_peak = max(max(abs(angle_range_map)));
    threshold_peak = max_peak * ratio;
    for rr = 1:size(angle_range_map, 1)
        if max(abs(angle_range_map(rr, :))) < threshold_peak
            continue
        end
        [peaks, aa] = findpeaks(abs(angle_range_map(rr, :)));
        for pp = 1:length(peaks)
            if aa(pp) == 1 || aa(pp) == size(angle_range_map, 2)
                continue;
            end
            if peaks(pp) >= threshold_peak
                angle_detect = [angle_detect; rr, aa(pp)];
            end
        end
    end
end
