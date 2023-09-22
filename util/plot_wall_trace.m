function plot_wall_trace(start, end_, angle)
    angle_matrix = [cos(angle) -sin(angle); 
        sin(angle) cos(angle)];
    start = start * angle_matrix';
    end_  = end_ * angle_matrix';
    % »æÖÆÂ·¾¶
    plot([start(1), end_(1)], [start(2), end_(2)], 'LineWidth',5)
end