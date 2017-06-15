close all
clear all
clc

% Constants
T = 0.01;
g = 9.81;
gyro_gauss = 0.001;  % Check datasheets for precise values
arw = 0.0015;
bias = 0.02;
acc_gauss = 0.3;
acc_rand = 0.1;
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
w_n = w + gyro_gauss*randn(1,length(w)) + arw*sqrt(1:length(w)) + bias;
ax_n = ax + acc_gauss*randn(1,length(ax)) + acc_rand*rand(1,length(ax).*cos(theta));
%az_n = az + acc_gauss*randn(1,length(az)) + acc_rand*rand(1,length(az).*sin(theta));

% Estimated thetas and differences
theta_est_gyro = zeros(1,length(w));
for i = 1:length(w)-1
    theta_est_gyro(i+1) = theta_est_gyro(i)+T*w_n(i);
end
theta_diff_gyro = (theta-theta_est_gyro);

theta_est_acc = asin(ax_n/g);
theta_diff_acc = theta-theta_est_acc;


figure()
subplot(2,2,1), plot(theta_est_gyro), title('Estimate based on gyroscope alone'), xlabel('time'), ylabel('rad');
subplot(2,2,2), plot(theta_diff_gyro), title('Error of gyroscope'), xlabel('time'), ylabel('rad');
subplot(2,2,3), plot(theta_est_acc), title('Estimate based on accelerometer alone'), xlabel('time'), ylabel('rad');
subplot(2,2,4), plot(theta_diff_acc), title('Error of accelerometer'), xlabel('time'), ylabel('rad');

value = [];
index = [];

%% Kalman filter 1:1

gyro_acc_sampling_ratio = 1;

% Without bias estimation, P modified every measurement 
theta_curr = 0;
theta_est_kalman1 = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        theta_curr = theta_curr + T*w_n((i-1)*gyro_acc_sampling_ratio+j);
        P = P + Q;
        theta_est_kalman1((i-1)*gyro_acc_sampling_ratio+j+1) = theta_curr;
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-theta_curr;
    S = P+R;
    G = P/S;
    theta_curr = theta_curr + G*v;
    P = (eye(size(G*g,1))-G*g)*P;
    
end

% Without bias estimation, P modified only during update phase
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

% With bias estimation, P modified every measurement
x =[0; 0];
A = [1 T;0 1];
B = [T; 0];
C = [1 0];
theta_est_kalman_bias1 = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        x = A*x + B*w_n((i-1)*gyro_acc_sampling_ratio+j);
        theta_est_kalman_bias1((i-1)*gyro_acc_sampling_ratio+j+1) = x(1);
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-C*x;
    S = C*P*C'+R;
    G = A*P*C'/S;
    x = x + G*v;
    P = A*P*A'-G*C*P*A'+Q;
    
end

% With bias estimation, P modified only during update phase
x =[0; 0];
A = [1 T;0 1];
B = [T; 0];
C = [1 0];
theta_est_kalman_bias2 = zeros(1,length(w)+1);
P = gyro_gauss;
Q = gyro_gauss;
R = acc_gauss;

for i = 1:(length(ax_n)/gyro_acc_sampling_ratio)
    % Prediction
    for j = 1:gyro_acc_sampling_ratio
        x = A*x + B*w_n((i-1)*gyro_acc_sampling_ratio+j);
        theta_est_kalman_bias2((i-1)*gyro_acc_sampling_ratio+j+1) = x(1);
    end
    
    % Update
    v = ax_n(i*gyro_acc_sampling_ratio)/g-C*x;
    S = C*P*C'+R;
    G = A*P*A'*C'/S;
    x = x + G*v;
    P = A*P*A'-G*C*A*P*A'+Q;
    
end

% Figures
figure()
subplot(2,1,1)
hold on
plot(theta), plot(theta_est_kalman1, 'r--'), plot(theta_est_kalman2, 'g.')
title('Real and Kalman estimated theta values without bias estimation'), xlabel('time'), ylabel('theta(rad)')
legend('real value','P every meas.', 'P update only')
hold off
subplot(2,1,2)
hold on
plot(theta), plot(theta_est_kalman_bias1, 'r--'), plot(theta_est_kalman_bias2, 'g.')
title('Real and Kalman estimated theta values with bias estimation'), xlabel('time'), ylabel('theta(rad)')
legend('real value','P every meas.', 'P update only')
hold off

err1 = theta-theta_est_kalman1(1:end-1);
err2 = theta-theta_est_kalman2(1:end-1);
err3 = theta-theta_est_kalman_bias1(1:end-1);
err4 = theta-theta_est_kalman_bias2(1:end-1);

figure()
subplot(2,2,1), plot(err1)
subplot(2,2,2), plot(err2)
subplot(2,2,3), plot(err3)
subplot(2,2,4), plot(err4)

[v, i] = min([err1*err1',err2*err2',err3*err3',err4*err4']);
value = [value, v];
index = [index, i];

