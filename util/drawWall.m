% ������ֻ���ڻ���שǽ
% ���룺����שǽ��ʼ�㡢����שǽ�����㡢����שǽ���
% �������
function drawWall(WallStartLoc, ...   % 1 x 2
    WallEndLoc, ...                   % 1 x 2
    WallThickness)                    % single number
    startLocation = [WallStartLoc(1) WallEndLoc(1)]; % ǽ�ڵ���ʼλ��
    endLocation   = [WallStartLoc(2) WallEndLoc(2)]; % ǽ�ڵĽ���λ��
    % �������ǽ��
    RelatedLocation = WallEndLoc - WallStartLoc;     % ǽ�ڵ����λ��
    wallAngle       = atan2(RelatedLocation(2), RelatedLocation(1)); % ǽ�ڵĽǶ�
    
    plot(startLocation, endLocation, 'LineWidth', WallThickness);
end