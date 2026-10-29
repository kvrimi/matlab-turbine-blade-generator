%% Turbine blade section modelling (original algorithm)
clear; clc; close all;

%% 1. 输入参数设置
r1 = 1;            % 叶片前缘半径 (mm)
r2 = 0.8;            % 叶片后缘半径 (mm)
rmax = 1.5;          % 叶片最大厚度半径 (mm)
alfa = 0*pi/180;  % 进口角 (rad)
beta = 55*pi/180;  % 出口角 (rad)
R = 7;           % 圆弧半径 (mm)

%% 2. 计算关键点坐标
% 点A: 原点
Ax = 0;
Ay = 0;

% 点B: 第一段直线终点
Bx = -3*tan(alfa);
By = -3;

% 点O: 圆弧中心
syms x y
eq1 = (x-Bx)^2 + (y-By)^2 == R^2;
eq2 = y-By == -cot(pi/2-alfa)*(x-Bx);
sol = solve([eq1, eq2], [x, y]);
x_sols = double(sol.x);
y_sols = double(sol.y);

if y_sols(1) > tan(alfa)*x_sols(1)
    Ox = x_sols(1);
    Oy = y_sols(1);
else
    Ox = x_sols(2);
    Oy = y_sols(2);
end

% 点C: 圆弧段终点
Cx = Ox + R*cos(beta);
Cy = Oy - R*sin(beta);

% 点D: 出口直线段终点
Dx = Cx + tan(beta)*((-11.2)-Cy);
Dy = -11.2;

fprintf('关键点坐标:\n');
fprintf('  A: (%.4f, %.4f)\n', Ax, Ay);
fprintf('  B: (%.4f, %.4f)\n', Bx, By);
fprintf('  O: (%.4f, %.4f)\n', Ox, Oy);
fprintf('  C: (%.4f, %.4f)\n', Cx, Cy);
fprintf('  D: (%.4f, %.4f)\n', Dx, Dy);

%% 3. 计算中弧线各段（高密度采样确保平滑）
% 3.1 AB段: 直线段（从A到B，不包括端点以便后续连接）
n_AB = 100;
y_AB = linspace(0, By, n_AB);
x_AB = y_AB * tan(alfa);

% 3.2 BC段: 圆弧段
angle_B = atan2(By-Oy, Bx-Ox);
angle_C = atan2(Cy-Oy, Cx-Ox);
n_BC = 200;
theta_BC = linspace(angle_B, angle_C, n_BC);
x_BC = Ox + R*cos(theta_BC);
y_BC = Oy + R*sin(theta_BC);

% 3.3 CD段: 直线段
n_CD = 100;
y_CD = linspace(Cy, Dy, n_CD);
x_CD = Cx + tan(beta)*(y_CD - Cy);

% 合并中弧线（去除重复点）
x_mid = [x_AB, x_BC(2:end), x_CD(2:end)];
y_mid = [y_AB, y_BC(2:end), y_CD(2:end)];

%% 4. 计算叶片厚度分布（平滑过渡）
n_points = length(y_mid);
r_thick = zeros(1, n_points);

% 使用三次样条插值实现平滑的厚度过渡
idx_max = round(n_points/3);
t_control = [1, idx_max, n_points];  % 控制点位置
r_control = [r1, rmax, r2];          % 控制点厚度
t_all = 1:n_points;
r_thick = interp1(t_control, r_control, t_all, 'pchip');  % 分段三次Hermite插值

%% 5. 生成叶盆曲线(Cp)和叶背曲线(Cq) - 真实包络线算法(完美相切)
x_Cp = zeros(1, n_points);
y_Cp = zeros(1, n_points);
x_Cq = zeros(1, n_points);
y_Cq = zeros(1, n_points);

for i = 1:n_points
    xn = x_mid(i);
    yn = y_mid(i);
    r = r_thick(i);

    % 1. 计算切向量和厚度差 (使用中心差分提高精度)
    if i == 1
        dx = x_mid(i+1) - x_mid(i);
        dy = y_mid(i+1) - y_mid(i);
        dr = r_thick(i+1) - r_thick(i);
    elseif i == n_points
        dx = x_mid(i) - x_mid(i-1);
        dy = y_mid(i) - y_mid(i-1);
        dr = r_thick(i) - r_thick(i-1);
    else
        dx = x_mid(i+1) - x_mid(i-1);
        dy = y_mid(i+1) - y_mid(i-1);
        dr = r_thick(i+1) - r_thick(i-1);
    end

    % 2. 归一化切向量 (T) 和法向量 (N)
    ds = sqrt(dx^2 + dy^2);
    tx = dx / ds;
    ty = dy / ds;

    % 法向量（逆时针旋转90度）
    nx = -ty;
    ny = tx;

    % 3. 计算厚度变化率 dr/ds
    dr_ds = dr / ds;
    % 限制范围防范数值截断误差导致开方出现复数
    dr_ds = max(min(dr_ds, 0.999), -0.999);

    % 4. 计算真实的包络线触点方向向量 (核心修正)
    % 数学推导：包络线方向向量 U = -dr/ds * T ± sqrt(1 - (dr/ds)^2) * N

    % 叶背(Cq)方向: 对应正法向 (+)
    ux_Cq = -dr_ds * tx + sqrt(1 - dr_ds^2) * nx;
    uy_Cq = -dr_ds * ty + sqrt(1 - dr_ds^2) * ny;

    % 叶盆(Cp)方向: 对应负法向 (-)
    ux_Cp = -dr_ds * tx - sqrt(1 - dr_ds^2) * nx;
    uy_Cp = -dr_ds * ty - sqrt(1 - dr_ds^2) * ny;

    % 5. 沿包络线方向生成表面坐标
    x_Cp(i) = xn + r * ux_Cp;
    y_Cp(i) = yn + r * uy_Cp;

    x_Cq(i) = xn + r * ux_Cq;
    y_Cq(i) = yn + r * uy_Cq;
end
%% 6. 生成前缘曲线(Ca)和后缘曲线(Cb) - 精确连接
% 前缘圆弧：从叶背起点到叶盆起点
angle_Cp_start = atan2(y_Cp(1) - Ay, x_Cp(1) - Ax);
angle_Cq_start = atan2(y_Cq(1) - Ay, x_Cq(1) - Ax);

% 确保角度范围正确（逆时针从叶背到叶盆）
if angle_Cp_start < angle_Cq_start
    angle_Cq_start = angle_Cq_start - 2*pi;
end

n_Ca = 100;
theta_Ca = linspace(angle_Cq_start, angle_Cp_start, n_Ca);
x_Ca = Ax + r1 * cos(theta_Ca);
y_Ca = Ay + r1 * sin(theta_Ca);

% 后缘圆弧：从叶盆终点到叶背终点
angle_Cp_end = atan2(y_Cp(end) - Dy, x_Cp(end) - Dx);
angle_Cq_end = atan2(y_Cq(end) - Dy, x_Cq(end) - Dx);

% 确保角度范围正确（逆时针从叶盆到叶背）
if angle_Cq_end < angle_Cp_end
    angle_Cq_end = angle_Cq_end + 2*pi;
end

n_Cb = 100;
theta_Cb = linspace(angle_Cp_end, angle_Cq_end, n_Cb);
x_Cb = Dx + r2 * cos(theta_Cb);
y_Cb = Dy + r2 * sin(theta_Cb);

%% 7. 组合完整轮廓（确保闭合）
% 逆时针顺序：Ca -> Cp -> Cb -> Cq(反向)
x_profile = [x_Ca, x_Cp, x_Cb, fliplr(x_Cq), x_Ca(1)];  % 闭合曲线
y_profile = [y_Ca, y_Cp, y_Cb, fliplr(y_Cq), y_Ca(1)];

%% 8. 绘图显示
figure('Position', [100, 100, 800, 800]);
hold on; grid on; axis equal;

% 主图：完整轮廓
subplot(1,1,1);
hold on; grid on; axis equal;
plot(x_profile, y_profile, 'b-', 'LineWidth', 2.5);
plot(x_mid, y_mid, 'r--', 'LineWidth', 1.5);
plot([Ax, Bx, Cx, Dx], [Ay, By, Cy, Dy], 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(Ax+1, Ay+1, 'A', 'FontSize', 12, 'FontWeight', 'bold');
text(Bx+1, By+1, 'B', 'FontSize', 12, 'FontWeight', 'bold');
text(Cx+1, Cy+1, 'C', 'FontSize', 12, 'FontWeight', 'bold');
text(Dx+1, Dy+1, 'D', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('X (mm)', 'FontSize', 12);
ylabel('Y (mm)', 'FontSize', 12);
title('涡轮叶片截面轮廓', 'FontSize', 14, 'FontWeight', 'bold');
legend('叶片轮廓', '中弧线', '关键点', 'Location', 'best');


%% 9. 曲线质量检查
fprintf('\n=== 曲线连续性检查 ===\n');
% 检查前缘连接
dist_Ca_Cp = sqrt((x_Ca(end)-x_Cp(1))^2 + (y_Ca(end)-y_Cp(1))^2);
fprintf('前缘到叶盆距离: %.8f mm\n', dist_Ca_Cp);

% 检查后缘连接
dist_Cp_Cb = sqrt((x_Cp(end)-x_Cb(1))^2 + (y_Cp(end)-y_Cb(1))^2);
fprintf('叶盆到后缘距离: %.8f mm\n', dist_Cp_Cb);

dist_Cb_Cq = sqrt((x_Cb(end)-x_Cq(end))^2 + (y_Cb(end)-y_Cq(end))^2);
fprintf('后缘到叶背距离: %.8f mm\n', dist_Cb_Cq);

dist_Cq_Ca = sqrt((x_Cq(1)-x_Ca(1))^2 + (y_Cq(1)-y_Ca(1))^2);
fprintf('叶背到前缘距离: %.8f mm\n', dist_Cq_Ca);

% 检查曲率连续性（通过相邻点的角度变化）
dx_profile = diff(x_profile);
dy_profile = diff(y_profile);
angles = atan2(dy_profile, dx_profile);
angle_changes = abs(diff(angles));
angle_changes(angle_changes > pi) = 2*pi - angle_changes(angle_changes > pi);
max_angle_change = max(angle_changes) * 180/pi;
fprintf('最大角度变化: %.4f 度\n', max_angle_change);
%% 9.5 简化点数（均匀采样）
fprintf('\n=== 减少点数 ===\n');
fprintf('原始点数: %d\n', length(x_profile));

% 设置目标点数
target_points = 500;  % 可以调整这个数字

% 均匀采样
indices = round(linspace(1, length(x_profile), target_points));
x_profile_simple = x_profile(indices);
y_profile_simple = y_profile(indices);

fprintf('简化后点数: %d\n', length(x_profile_simple));
fprintf('压缩比: %.1f%%\n', (1-length(x_profile_simple)/length(x_profile))*100);

% 验证数据
fprintf('\n=== 数据验证 ===\n');
fprintf('x_profile_simple 长度: %d\n', length(x_profile_simple));
fprintf('y_profile_simple 长度: %d\n', length(y_profile_simple));

%% 10. 创建Excel表格（使用简化数据）
fprintf('\n=== 创建Excel表格 ===\n');

% 确保使用简化后的数据
T = table();
T.X = x_profile_simple(:);  % 注意这里用的是 x_profile_simple
T.Y = y_profile_simple(:);  % 注意这里用的是 y_profile_simple
T.Z = zeros(length(x_profile_simple), 1);

% 验证表格行数
fprintf('表格行数: %d\n', height(T));

% 保存为Excel
filename = 'blade_profile_simple.xlsx';
writetable(T, filename);

fprintf('✓ Excel文件已生成: %s\n', filename);
fprintf('✓ 表格行数: %d 行\n', height(T));
