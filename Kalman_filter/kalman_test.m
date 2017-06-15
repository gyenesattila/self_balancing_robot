close all
clear all
clc

% Constants
T = 0.01;
gyro_acc_sampling_ratio = 1;
g = 9.81;
gyro_gauss = 0.001;  % Check datasheets for precise values
arw = 0;
acc_gauss = 0.3;
acc_rand = 0;
w_sensitivity = 100*(pi/180);

% Exact values
w = w_sensitivity*rand(1,1000)-w_sensitivity/2;
theta = zeros(1,length(w));

for i = 1:length(w)-1
    theta(i+1) = theta(i)+T*w(i);
end

ax = sin(theta).*g;
az = cos(theta).*g;

% Added noise
w_n = w + gyro_gauss*randn(1,length(w)) + arw*sqrt(1:length(w))+0,02;
ax_n = ax + acc_gauss*randn(1,length(ax)) + acc_rand*rand(1,length(ax).*cos(theta));
az_n = az + acc_gauss*randn(1,length(az))% + acc_rand*rand(1,length(az).*sin(theta));

% Estimated thetas and differences
theta_est_gyro = zeros(1,length(w));
for i = 1:length(w)-1
    theta_est_gyro(i+1) = theta_est_gyro(i)+T*w_n(i);
end
theta_diff_gyro = (theta-theta_est_gyro);
theta_diff_gyro_percent = 100*(theta-theta_est_gyro)./theta;

theta_est_acc = asin(ax_n/g);
theta_diff_acc = theta-theta_est_acc;
theta_diff_acc_percent = 100*(theta-theta_est_acc)./theta;

% a_sum = ax.*az;
% theta_est_acc2 = asin(ax_n/g);
% theta_diff_acc2 = theta-theta_est_acc2;
% theta_diff_acc2_percent = 100*(theta-theta_est_acc2)./theta;

figure(1)
subplot(2,3,1), plot(theta_est_gyro);
subplot(2,3,2), plot(theta_diff_gyro);
subplot(2,3,3), plot(theta_diff_gyro_percent);
subplot(2,3,4), plot(theta_est_acc);
subplot(2,3,5), plot(theta_diff_acc);
subplot(2,3,6), plot(theta_diff_acc_percent);

%% Kalman filter

theta_curr = 0;
theta_est_kalman = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        theta_curr = theta_curr + T*w_n((i-1)*gyro_acc_sampling_ratio+j);
        P = P + Q;
        theta_est_kalman((i-1)*gyro_acc_sampling_ratio+j+1) = theta_curr;
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-theta_curr;
    S = P+R;
    G = P/S;
    theta_curr = theta_curr + G*v;
    P = (eye(size(G*g,1))-G*g)*P;
    
end

theta_curr = 0;
theta_est_kalman2 = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        theta_curr = theta_curr + T*w_n((i-1)*gyro_acc_sampling_ratio+j);
        %P = P + Q;
        theta_est_kalman2((i-1)*gyro_acc_sampling_ratio+j+1) = theta_curr;
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-theta_curr;
    S = P+R;
    G = P/S;
    theta_curr = theta_curr + G*v;
    P = (eye(size(G*g,1))-G*g)*P+Q;
    
end

figure(2)
hold on
plot(theta)
plot(theta_est_kalman, 'r--')
%plot(theta_est_kalman2, 'g.')
hold off

figure(3)
subplot(2,2,1), plot(theta-theta_est_gyro)
subplot(2,2,2), plot(theta-theta_est_acc)
subplot(2,2,3), plot(theta-theta_est_kalman(1:end-1))
subplot(2,2,4), plot(theta-theta_est_kalman2(1:end-1))

% figure(4)
% plot(theta_est_kalman-theta_est_kalman2)

x =[0; 0];
A = [1 T;0 1];
B = [T; 0];
C = [1 0];
theta_est_kalman = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        x = A*x + B*w_n((i-1)*gyro_acc_sampling_ratio+j);
        theta_est_kalman((i-1)*gyro_acc_sampling_ratio+j+1) = x(1);
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-C*x;
    S = C*P*C'+R;
    G = A*P*C'/S;
    x = x + G*v;
    P = A*P*A'-G*C*P*A'+Q;
    
end

figure(5)
subplot(1,2,1), plot(theta-theta_est_kalman(1:end-1))
subplot(1,2,2), plot(theta-theta_est_kalman2(1:end-1))


x =[0; 0];
A = [1 T;0 1];
B = [T; 0];
C = [1 0];
theta_est_kalman2 = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        x = A*x + B*w_n((i-1)*gyro_acc_sampling_ratio+j);
        theta_est_kalman2((i-1)*gyro_acc_sampling_ratio+j+1) = x(1);
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-C*x;
    S = C*P*C'+R;
    G = A*P*A'*C'/S;
    x = x + G*v;
    P = A*P*A'-G*C*A*P*A'+Q;
    
end

figure(6)
subplot(1,3,1), plot(theta-theta_est_kalman(1:end-1))
subplot(1,3,2), plot(theta-theta_est_kalman2(1:end-1))
subplot(1,3,3), plot(theta_est_kalman-theta_est_kalman2)
