function plot_imu_data(imu_data, t)
% PLOT_IMU_DATA 可视化 6 轴 IMU 数据 Deepseek 生成可视化代码
%   plot_imu_data(imu_data, t)
%
% 输入:
%   imu_data - 6×N 矩阵，前 3 行为陀螺角增量 (x,y,z)，后 3 行为加计速度增量 (x,y,z)
%   t        - 时间信息，可以是：
%              * 1×N 时间向量（单位：秒）
%              * 标量采样间隔 dt（将自动生成 0:dt:(N-1)*dt）
%              若省略，则默认 t = 1:N
%
% 输出:
%   无返回值，直接绘制图形

    % ---------- 输入参数检查 ----------
    if nargin < 1
        error('至少需要提供 IMU 数据矩阵。');
    end
    
    [nRows, nCols] = size(imu_data);
    if nRows ~= 6
        error('IMU 数据必须为 6×N 矩阵（前 3 行陀螺，后 3 行加计）。');
    end
    
    % 处理时间向量
    if nargin < 2 || isempty(t)
        t = 1:nCols;
    elseif isscalar(t)
        dt = t;
        t = (0:nCols-1) * dt;
    else
        t = t(:)';   % 确保为行向量
        if length(t) ~= nCols
            error('时间向量长度必须与数据列数一致。');
        end
    end
    
    % ---------- 提取各轴数据 ----------
    gyro_x = imu_data(1, :);
    gyro_y = imu_data(2, :);
    gyro_z = imu_data(3, :);
    acc_x  = imu_data(4, :);
    acc_y  = imu_data(5, :);
    acc_z  = imu_data(6, :);
    
    % ---------- 创建图形窗口 ----------
    figure('Name', 'IMU Data Visualization', 'Color', 'w', ...
           'Units', 'centimeters', 'Position', [2, 2, 24, 12]);
    
    % 定义颜色
    colors = [0.00 0.45 0.74;   % 蓝
              0.85 0.33 0.10;   % 橙
              0.47 0.67 0.19];  % 绿
    lineWidth = 1.2;
    
    % ---------- 第一行：陀螺仪 ----------
    subplot(2, 3, 1);
    plot(t, gyro_x, 'Color', colors(1,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Angular Increment (rad)');
    title('Gyro X'); grid on; set(gca, 'FontSize', 10);
    
    subplot(2, 3, 2);
    plot(t, gyro_y, 'Color', colors(2,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Angular Increment (rad)');
    title('Gyro Y'); grid on; set(gca, 'FontSize', 10);
    
    subplot(2, 3, 3);
    plot(t, gyro_z, 'Color', colors(3,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Angular Increment (rad)');
    title('Gyro Z'); grid on; set(gca, 'FontSize', 10);
    
    % ---------- 第二行：加速度计 ----------
    subplot(2, 3, 4);
    plot(t, acc_x, 'Color', colors(1,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Velocity Increment (m/s)');
    title('Accel X'); grid on; set(gca, 'FontSize', 10);
    
    subplot(2, 3, 5);
    plot(t, acc_y, 'Color', colors(2,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Velocity Increment (m/s)');
    title('Accel Y'); grid on; set(gca, 'FontSize', 10);
    
    subplot(2, 3, 6);
    plot(t, acc_z, 'Color', colors(3,:), 'LineWidth', lineWidth);
    xlabel('Time (s)'); ylabel('Velocity Increment (m/s)');
    title('Accel Z'); grid on; set(gca, 'FontSize', 10);
    
    % 总标题
    sgtitle('6-Axis IMU Data Visualization when align', 'FontSize', 14, 'FontWeight', 'bold');
end