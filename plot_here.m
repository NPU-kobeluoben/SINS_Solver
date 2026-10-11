function  plot_here(avp,avp_buff_results)
%% 捷联惯导解算结果与参考 AVP 对比可视化 (Deepseek AI生成可视化代码)
t_ref = avp(:,10);
t_sins = t_ref;   % 假设两者时间戳一致；若不一致，请调整

% 参考值提取
att_ref  = avp(:,1:3);    % [pitch roll yaw] (rad)
vel_ref  = avp(:,4:6);    % [vE vN vU] (m/s)
pos_ref  = avp(:,7:9);    % [lat lon h] (rad, rad, m)

% SINS 结果提取（转置为 N×9）
att_sins = avp_buff_results(1:3,:)';
vel_sins = avp_buff_results(4:6,:)';
pos_sins = avp_buff_results(7:9,:)';

% 角度转为度（便于观察）
att_ref_deg  = rad2deg(att_ref);
att_sins_deg = rad2deg(att_sins);

% 位置纬经转为度
pos_ref_deg  = pos_ref;  pos_ref_deg(:,1:2)  = rad2deg(pos_ref(:,1:2));
pos_sins_deg = pos_sins; pos_sins_deg(:,1:2) = rad2deg(pos_sins(:,1:2));

% 处理方位角跳变（unwrap，避免 0/360 处跳变）
att_ref_deg(:,3)  = unwrap(att_ref_deg(:,3)  * pi/180) * 180/pi;
att_sins_deg(:,3) = unwrap(att_sins_deg(:,3) * pi/180) * 180/pi;

% ======================== 图像1：状态对比 ========================
figure('Name','SINS vs Reference AVP','Color','w','Position',[100 50 1200 800]);

% ---- (1) 俯仰角 ----
subplot(3,3,1);
plot(t_ref, att_ref_deg(:,1), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, att_sins_deg(:,1), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Pitch (deg)');
title('俯仰角对比'); legend('Ref','SINS','Location','best');

% ---- (2) 横滚角 ----
subplot(3,3,2);
plot(t_ref, att_ref_deg(:,2), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, att_sins_deg(:,2), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Roll (deg)');
title('横滚角对比'); legend('Ref','SINS','Location','best');

% ---- (3) 方位角 ----
subplot(3,3,3);
plot(t_ref, att_ref_deg(:,3), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, att_sins_deg(:,3), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Yaw (deg)');
title('方位角对比'); legend('Ref','SINS','Location','best');

% ---- (4) 东向速度 ----
subplot(3,3,4);
plot(t_ref, vel_ref(:,1), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, vel_sins(:,1), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('v_E (m/s)');
title('东向速度对比'); legend('Ref','SINS','Location','best');

% ---- (5) 北向速度 ----
subplot(3,3,5);
plot(t_ref, vel_ref(:,2), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, vel_sins(:,2), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('v_N (m/s)');
title('北向速度对比'); legend('Ref','SINS','Location','best');

% ---- (6) 天向速度 ----
subplot(3,3,6);
plot(t_ref, vel_ref(:,3), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, vel_sins(:,3), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('v_U (m/s)');
title('天向速度对比'); legend('Ref','SINS','Location','best');

% ---- (7) 纬度 ----
subplot(3,3,7);
plot(t_ref, pos_ref_deg(:,1), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, pos_sins_deg(:,1), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Latitude (deg)');
title('纬度对比'); legend('Ref','SINS','Location','best');

% ---- (8) 经度 ----
subplot(3,3,8);
plot(t_ref, pos_ref_deg(:,2), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, pos_sins_deg(:,2), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Longitude (deg)');
title('经度对比'); legend('Ref','SINS','Location','best');

% ---- (9) 高度 ----
subplot(3,3,9);
plot(t_ref, pos_ref(:,3), 'k-', 'LineWidth', 1.5); hold on;
plot(t_sins, pos_sins(:,3), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Height (m)');
title('高度对比'); legend('Ref','SINS','Location','best');

sgtitle('SINS 解算结果 vs 参考 AVP','FontSize', 14, 'FontWeight', 'bold');
%% ======================== 图像2：误差曲线 ========================
figure('Name','SINS Errors','Color','w','Position',[150 80 1200 800]);

% 误差 = SINS - Ref
err_att = att_sins_deg - att_ref_deg;                 % deg
err_vel = vel_sins - vel_ref;                         % m/s
err_pos = pos_sins_deg - pos_ref_deg;                 % deg, deg, m

% ---- (1) 俯仰角误差 ----
subplot(3,3,1);
plot(t_sins, err_att(:,1), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta Pitch (deg)');
title('俯仰角误差');

% ---- (2) 横滚角误差 ----
subplot(3,3,2);
plot(t_sins, err_att(:,2), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta Roll (deg)');
title('横滚角误差');

% ---- (3) 方位角误差 ----
subplot(3,3,3);
plot(t_sins, err_att(:,3), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta Yaw (deg)');
title('方位角误差');

% ---- (4) 东向速度误差 ----
subplot(3,3,4);
plot(t_sins, err_vel(:,1), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta v_E (m/s)');
title('东向速度误差');

% ---- (5) 北向速度误差 ----
subplot(3,3,5);
plot(t_sins, err_vel(:,2), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta v_N (m/s)');
title('北向速度误差');

% ---- (6) 天向速度误差 ----
subplot(3,3,6);
plot(t_sins, err_vel(:,3), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta v_U (m/s)');
title('天向速度误差');

% ---- (7) 纬度误差 ----
subplot(3,3,7);
plot(t_sins, err_pos(:,1), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta Lat (deg)');
title('纬度误差');

% ---- (8) 经度误差 ----
subplot(3,3,8);
plot(t_sins, err_pos(:,2), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta Lon (deg)');
title('经度误差');

% ---- (9) 高度误差 ----
subplot(3,3,9);
plot(t_sins, err_pos(:,3), 'b-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('\Delta h (m)');
title('高度误差');

sgtitle('SINS 解算误差（SINS - Reference）','FontSize', 14, 'FontWeight', 'bold');

end

