function drawUOIarea(ra_map, detect_result, range_axis, angle_axis)
% detect_result: 1: rangle 2:angle
    
figure(11110)
imagesc(angle_axis, range_axis, ra_map)
hold on
xlabel("Degree (rad)")
ylabel("Distance (meter)")
scatter(detect_result(2,:), detect_result(1,:), 'r')
legend('Detected UOI')

end

