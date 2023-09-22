% 本函数只用于绘制砖墙
% 输入：反射砖墙起始点、反射砖墙结束点、反射砖墙厚度
% 输出：无
function drawWall(WallStartLoc, ...   % 1 x 2
    WallEndLoc, ...                   % 1 x 2
    WallThickness)                    % single number
    startLocation = [WallStartLoc(1) WallEndLoc(1)]; % 墙壁的起始位置
    endLocation   = [WallStartLoc(2) WallEndLoc(2)]; % 墙壁的结束位置
    % 方便绘制墙壁
    RelatedLocation = WallEndLoc - WallStartLoc;     % 墙壁的相对位置
    wallAngle       = atan2(RelatedLocation(2), RelatedLocation(1)); % 墙壁的角度
    
    plot(startLocation, endLocation, 'LineWidth', WallThickness);
end