%% Considering Fill-Factor (FF) for MIMO OFDM
clc
clear
close all 
% Define constants and parameters
z = 2;
d = z;
lambda = 850e-9; 
Pt = 1e-3; 
P_elec = (1/9)*Pt^2;
B = 20e9; % VCSEL BW 
RIN_db = -155; % dB/Hz, Laser noise 
RIN = 10^(RIN_db/10);
r_PD = 3e-3; % PD radius (single PD)
R_PD = 0.4; % responsivity 
A_PD = 28*1e-4; % PD area (single PD) 
RL = 50; % Load resistance 
Fn_db = 5; % TIA noise figure 
Fn = 10^(5/10); % Fn 
K_boltz = 1.38*1e-23; % Boltzmann 
T = 298; % Temperature in Kelvin 
q = 1.60217662*1e-19; % elementary charge (C) 
BER = 1e-3; % FEC limit 
Gamma = - log(5 * BER)/1.5;
% log2(M) % Number of transmitted bits per channel use 
% In the case of using OFDM: 
N = 256; % number of symbols N_FFT
symbol_rate = 2*B / N;
zeta = (N - 2)/N; % subcarrier utilization factor 
w_0 = 100e-6;
zR = (pi .* (w_0.^2)) ./ lambda; 
wL = w_0 .* sqrt(1 + (z ./ zR).^2); % Beam waist at distance z

Nt = 25; 
% Nr3 = 81; % Nr (PDs) for config III
Nr = [81];
r_PD = 3e-3; 
delta = 0e-3;
d_PD = 2*r_PD + delta; 

x_DE = (0:0.1:60).*1e-3; 
y_DE = zeros(1,length(x_DE));
r_DE = x_DE; 
syms x y; 
sqrt_pi = sqrt(pi);
sqrt_2 = sqrt(2);

figure
hold on 
ylabel('Rate [Tb/s]')
xlabel('r_D_E [mm]');
grid on 
for nr = 1:length(Nr)

    Hij_MIMO_dis = zeros(Nr(nr),Nt);
    H_MIMO_dis = zeros(Nr(nr),Nt);
    Hii_MIMO_dis = zeros(Nr(nr),1); % Note: Only diagonal elements are needed
    var_dis = zeros(1, Nr(nr));
    SNR_dis = zeros(1, Nr(nr));
    Rate_dis = zeros(1, Nr(nr));

    for k = 1:length(r_DE)
        for i = 1:Nr(nr)
            K = ceil(sqrt(Nr(nr)));
            mi = floor((i-1)/K); % row in H matrix
            ni = i - mi * K; % column in H matrix
            xi = (- (K - 1)/2 + ni - 1)*d_PD;
            yi = ((K - 1)/2 - mi)*d_PD;

            % mi = floor(i/K); % row in H matrix
            % ni = i - (floor(i/K) - 1)*K; % column in H matrix
            % xi = (- (K - 1)/2 + ni - 1)*d_PD;
            % yi = ( (K - 1)/2 + mi + 1)*d_PD;

            for j = 1:Nt
                Nt_sq = ceil(sqrt(Nt));
                mj = floor((j-1)/Nt_sq); % row in H matrix
                nj = j - mj * Nt_sq; % column in H matrix
                xj = (- (Nt_sq - 1)/2 + nj - 1)*d_PD;
                yj = ((Nt_sq - 1)/2 - mj)*d_PD;
                
                % mj = floor(j/Nt_sq); % row in H matrix
                % nj = j - (floor(j/Nt_sq) - 1)*Nt_sq; % column in H matrix
                % xj = (- (Nt_sq - 1)/2 + nj - 1)*d_PD;
                % yj = ( (Nt_sq - 1)/2 + mj + 1)*d_PD;

                if i ~= j
                    arg_x1_dis = (sqrt_pi * r_PD + 2 * (xi - xj - x_DE(k))) / (sqrt_2 * wL);
                    arg_x2_dis = (sqrt_pi * r_PD - 2 * (xi - xj - x_DE(k))) / (sqrt_2 * wL);
                    arg_y1_dis = (sqrt_pi * r_PD + 2 * (yi - yj - y_DE(k))) / (sqrt_2 * wL);
                    arg_y2_dis = (sqrt_pi * r_PD - 2 * (yi - yj - y_DE(k))) / (sqrt_2 * wL);
                    term_x_dis = (erf(arg_x1_dis) + erf(arg_x2_dis));
                    term_y_dis = (erf(arg_y1_dis) + erf(arg_y2_dis));
                    Hij_MIMO_dis(i,j) = 0.25*term_x_dis* term_y_dis;
                    H_MIMO_dis(i,j) = Hij_MIMO_dis(i,j);
                else
                    arg_x1_dis = (sqrt_pi * r_PD + 2 * (- x_DE(k))) / (sqrt_2 * wL);
                    arg_x2_dis = (sqrt_pi * r_PD - 2 * (- x_DE(k))) / (sqrt_2 * wL);
                    arg_y1_dis = (sqrt_pi * r_PD + 2 * (- y_DE(k))) / (sqrt_2 * wL);
                    arg_y2_dis = (sqrt_pi * r_PD - 2 * (- y_DE(k))) / (sqrt_2 * wL);
                    term_x_dis = (erf(arg_x1_dis) + erf(arg_x2_dis));
                    term_y_dis = (erf(arg_y1_dis) + erf(arg_y2_dis));
                    Hii_MIMO_dis(i) = 0.25*term_x_dis* term_y_dis;
                    H_MIMO_dis(i,j) = Hii_MIMO_dis(i);
                end
            end

            var_dis(i) = (4*K_boltz*T/RL)*B*Fn + 2*q*(sum(R_PD*Hij_MIMO_dis(i,:)*Pt))*B + RIN*(sum((R_PD*Hij_MIMO_dis(i,:)*Pt).^2))*B;
            SNR_dis(i) = ((R_PD^2) * (Hii_MIMO_dis(i).^2) * P_elec) / (sum((R_PD^2).*(Hij_MIMO_dis(i,:).^2).*P_elec) + var_dis(i));
            Rate_dis(i) = zeta*B*log2(1 + SNR_dis(i)/Gamma);  % Transmission rate for VCSEL_i
        end

        [~,S_dis,~] = svd(H_MIMO_dis);
        S_dis_diag = diag(S_dis);
        threshold = 0.0001 * max(S_dis_diag);
        selected_diag_elements = S_dis_diag(S_dis_diag > threshold);
        lambdai = diag(diag(selected_diag_elements));
        var_dis_svd = (4*K_boltz*T/RL)*B*Fn + 2*q.*(R_PD* lambdai*Pt).*B + RIN.*((R_PD* lambdai*Pt).^2).*B;
        SNR_svd_dis = ((R_PD^2)*P_elec.*(lambdai.^2)) ./ var_dis_svd;
        Rate_svd_dis = zeta.*B.*log2(1 + SNR_svd_dis./Gamma);  % Transmission rate for VCSEL_j
        sum_Rate_dis(k) =  sum(Rate_dis(:))./1e12; % sum rate of MIMO OFDM (Tbps)
        sum_Rate_svd_dis(k) =  sum(Rate_svd_dis(:))./1e12; % sum rate of multiple identical SISO OFDM links (Tbps)

    end
    wosvd = sum_Rate_dis;
    wsvd = sum_Rate_svd_dis;
    plot(r_DE.*1e3,wosvd,'b','LineWidth',2)
    plot(r_DE.*1e3,wsvd,'r','LineWidth',2)
end