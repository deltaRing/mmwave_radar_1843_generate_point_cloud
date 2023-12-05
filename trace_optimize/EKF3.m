% 扩展卡尔曼滤波3 
% 输入1：雷达观测到的俯仰角 E
% 输入2：雷达观测到的方位角 A
% 输入3：雷达观测到的距离   R
% 输入4：过去航迹的观测值 Z (1 x 6) X Vx Y Vy Z Vz
% 输入5：过去航迹的协方差 P 
%  (过去的航迹观测大小和雷达新观测到的大小必须一致)
%  (如果没有观测到新的数据那就不要加进去（可以为NULL）)
% 可选输入1：EKF2 二维还是三维
% 可选输入2：T    单帧的
% 输出1：新的航迹观测值 Zn
% 输出2：新的航迹协方差 Pn
function [Zn, Pn] = EKF3(E, A, R, Zn, Pn, EKF2, t)
% 设置默认值
if nargin == 5, EKF2 = 0; t = 0.1; end
sz = -1;
% 如果是 EKF2 维度是 4 如果是 EKF3 维度是 6
if EKF2, sz = 2; sh = 2; else sz = 3; sh = 3; end
% 初始化噪声
Q_ = 10;     % 运动方程中的噪声
R_ = 1;      % 观测方程中的噪声
% 运动方程
% Xt = Xt-1 + Vxt-1 * t 
% Yt = Yt-1 + Vyt-1 * t
% Zt = Zt-1 + Vzt-1 * t
% Vx = Vxt-1 + n(t)
% Vy = Vyt-1 + n(t)
% Vz = Vzt-1 + n(t)
% 如何将运动方程转变为观测方程
% azi   = arctan(y/x)
% ele   = arctan(z/√(x^2 + y^2))
% range = √(x^2 + y^2 + z^2) 
% 运动方程求导
% [1 t 0 0 0 0;  X
%  0 1 0 0 0 0;  Vx
%  0 0 1 t 0 0;  Y
%  0 0 0 1 0 0;  Vy
%  0 0 0 0 1 t;  Z
%  0 0 0 0 0 1]; Vz
% 观测方程求导
% [ -y / (x^2 + y^2) 0 x / (x^2 + y^2) 0 0 0;
%    x * z / (range^2 * range_xy) 0 y * z / (range^2 * range_xy) 0 -range_xy / range^2 0;
%    x / range 0 y / range 0 z / range 0]
% 其中
% range_xy = √(x^2 + y^2)
% range = √(x^2 + y^2 + z^2)

for ii = 1:size(Zn, 1)
    if EKF2

    else
        % 如果有数据
        % 构造真实状态
        if ~isempty(E{ii})
            r_xy = R{ii} * cos(E{ii});
            r_z_  = R{ii} * sin(E{ii});
            r_x_  = r_xy  * cos(A{ii});
            r_y_  = r_xy  * sin(A{ii});

            r_x  = Zn(ii, 1) + Zn(ii, 2) * t;
            r_y  = Zn(ii, 3) + Zn(ii, 4) * t;
            r_z  = Zn(ii, 5) + Zn(ii, 6) * t;
            v_x  = (r_x_ - Zn(ii, 1)) / t;
            v_y  = (r_y_ - Zn(ii, 3)) / t;
            v_z  = (r_z_ - Zn(ii, 5)) / t;
            X_ = [r_x, v_x, r_y, v_y, r_z, v_z]; % 预测的状态
            Zf = [atan2(X_(3), X_(1)) ...
                    atan2(X_(5), norm([X_(1), X_(3)])) ...
                    norm([X_(1), X_(3), X_(5)])];
            Z_ = [A{ii}, E{ii}, R{ii}];          % 真实的观测值
        else
            % 不进行更新
            continue;
        end
        
        range = norm([r_x r_y r_z]);
        range_xy = norm([r_x r_y]);

        F = [1 t 0 0 0 0; 
            0 1 0 0 0 0;
            0 0 1 t 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 t;
            0 0 0 0 0 1]; % 6 x 6
        H = [-r_y / range_xy^2 0 r_x / range_xy^2 0 0 0;
            r_x * r_z / (range^2 * range_xy) 0 r_y * r_z / (range^2 * range_xy) 0 -range_xy / range^2 0;
            r_x / range 0 r_y / range 0 r_z / range 0]; % 6 x 3
    end
    
    P_ = F * Pn{ii} * F' + Q_ * eye(6); % 6 x 6
    K = P_ * H' / (H * P_ * H' + R_ * eye(3));
    X = X_ + (K * (Z_ - Zf)')';
    Pn{ii} = P_ - K * H * P_;          % 更新协方差
    Zn(ii, :) = X;                     % 更新状态
end
end