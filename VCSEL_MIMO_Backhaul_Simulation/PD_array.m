clc
clear
close all 
Nr = 25;
r_PD = 3e-3; 
delta = 6e-3; 
d_PD = 2 * r_PD + delta; 
K = ceil(sqrt(Nr)); 
for i = 1:Nr
    mi = floor((i - 1) / K);
    ni = i - mi * K;
    xi = (- (K - 1) / 2 + (ni - 1)) * d_PD;
    yi = ((K - 1) / 2 - mi) * d_PD;

    % mi = floor((i) / K); 
    % ni = i - (floor(i/K) - 1)*K;
    % xi = (- (K-1)/2 + ni - 1)*d_PD;
    % yi = ((K-1)/2 - mi + 1)*d_PD; 

    x_positions(i) = xi;
    y_positions(i) = yi;
end
x_positions = x_positions * 1000; 
y_positions = y_positions * 1000; 
r_PD_mm = r_PD * 1000;
figure;
hold on;
for i = 1:Nr
    rectangle('Position', [x_positions(i) - r_PD_mm, y_positions(i) - r_PD_mm, 2 * r_PD_mm, 2 * r_PD_mm], ...
              'Curvature', [1, 1], ...
              'FaceColor', 'g', ...
              'EdgeColor', 'k');
end
axis equal;
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
title('PD Array');
grid on;


Ky = K - 1; 
for i = 1:Ky^2
    mi = floor((i - 1) / Ky); 
    ni = i - mi * Ky; 
    xi = (- (Ky - 1) / 2 + (ni - 1)) * d_PD;
    yi = ((Ky - 1) / 2 - mi) * d_PD;
    x_y(i) = xi;
    y_y(i) = yi;
end
x_y = x_y*1000; 
y_y = y_y*1000; 
for i = 1:((K -1)^2)
    rectangle('Position', [x_y(i) - r_PD_mm, y_y(i) - r_PD_mm, 2 * r_PD_mm, 2 * r_PD_mm], ...
              'Curvature', [1, 1], ...
              'FaceColor', 'y', ...
              'EdgeColor', 'k');
end
