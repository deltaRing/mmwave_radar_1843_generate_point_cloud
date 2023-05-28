function angle_detect = AnglePeaksDetect(angle_range_map)
    angle_detect = [];
    range_index = [];
    if length(size(angle_range_map)) ~= 2
        fprintf("需要 距离-角度2维谱\n");
        return;
    end
    
    ratio = 0.35;
    max_peak = max(max(abs(angle_range_map)));
    threshold_peak = max_peak * ratio;
    r_gate_of_interest = std(angle_range_map');
    if (min(r_gate_of_interest) * 20.0 > max(r_gate_of_interest))
        fprintf("没有感兴趣的角度谱\n")
        return;
    end
    [value, index] = find(r_gate_of_interest > mean(r_gate_of_interest) * 0.5);
    for ii = 1:length(index)
        rr = index(ii);
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