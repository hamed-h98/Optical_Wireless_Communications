%% Displacement MIMO
% H0 (x_DE, y_DE)

% MIMO
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
Fn_db = 5; %TIA noise figure 
Fn = 10^(5/10); % Fn 
K_boltz = 1.38*1e-23; % Boltzmann 
T = 298; % Temperature in Kelvin 
q = 1.5723*1e-34; % elementary charger (C) 
BER = 1e-3; % FEC limit 
Gamma = - log(5 * BER)/1.5;
% log2(M) % Number of transmitted bits per channel use 
% In the case of using OFDM: 
N = 64; % number of symbols N_FFT
symbol_rate = 2*B / N;
zeta = (N - 2)/N; % subcarrier utilization factor 
w0 = 10e-6:1e-6:100e-6;  % Initial beam waist 
delta = 6e-3; % inter-element spacing in MIMO
% center-to-center distance for neighboring PDs along rows or columns of the array:
d_PD = 2*r_PD + delta; 
K = [4 9 16 25];
figure
hold on 
xlabel('\omega_0 [\mum]');
ylabel(['Rate [Tb/s]']);
grid on 

for l = 1:length(w0)
    zR(l) = (pi*(w0(l)^2))/lambda;
    wL_sq(l) = (w0(l).^2)*(1 + (z./zR(l))^2); % Beam waist squared at distance z (w(L)^2)
end

% Integral bounds: 
x_min = @(y) -sqrt(r_PD^2 - y.^2);
x_max = @(y) sqrt(r_PD^2 - y.^2);
y_min = - r_PD;
y_max = r_PD;

L = 2; % Dsitance between Tx and Rx 
Nt = 25; 
w_0 = 100e-6;
zR = (pi .* (w_0.^2)) ./ lambda; 
wL = w_0 .* sqrt(1 + (z ./ zR).^2); % Beam waist at distance z
% Nr = 25; % Nr (PDs) for config I, 20% fill factor (FF = Nr*A_PD/(W^2))
r_PD = 3e-3; 
delta = 6e-3;
% delta = 0;
d_PD = 2*r_PD + delta; 
Nr = 25; 
x_DE = (0:0.1:60).*1e-3; 
y_DE = zeros(1,length(x_DE));
r_DE = x_DE; 
syms x y; 
figure
hold on 
ylabel('Rate [Tb/s]')
xlabel('r_D_E [mm]');
grid on 
sqrt_pi = sqrt(pi);
sqrt_2 = sqrt(2);

for k = 1:length(r_DE)
    for i = 1:Nr
        K = sqrt(Nr);
        mi = floor(i/K); % row in H matrix
        ni = i - (floor(i/K) - 1)*K; % column in H matrix
        xi = (- (K - 1)/2 + ni - 1)*d_PD;
        yi = ( (K - 1)/2 + mi + 1)*d_PD;

        for j = 1:Nt
            Nr_sq = sqrt(Nt);
            mj = floor(j/Nr_sq); % row in H matrix
            nj = j - (floor(j/Nr_sq) - 1)*Nr_sq; % column in H matrix
            xj = (- (Nr_sq - 1)/2 + nj - 1)*d_PD;
            yj = ( (Nr_sq - 1)/2 + mj + 1)*d_PD;

            if i ~= j
                arg_x1_dis = (sqrt_pi * r_PD + 2 * (xi - xj - x_DE(k))) / (sqrt_2 * wL);
                arg_x2_dis = (sqrt_pi * r_PD - 2 * (xi - xj - x_DE(k))) / (sqrt_2 * wL);
                arg_y1_dis = (sqrt_pi * r_PD + 2 * (yi - yj - y_DE(k))) / (sqrt_2 * wL);
                arg_y2_dis = (sqrt_pi * r_PD - 2 * (yi - yj - y_DE(k))) / (sqrt_2 * wL);
                term_x_dis = (erf(arg_x1_dis) + erf(arg_x2_dis));
                term_y_dis = (erf(arg_y1_dis) + erf(arg_y2_dis));
                Hij_MIMO_dis(i,j) = 0.25*term_x_dis* term_y_dis;
                H_MIMO_dis(i,j) = Hij_MIMO_dis(i,j);

            else if i == j
                    arg_x1_dis = (sqrt_pi * r_PD + 2 * (- x_DE(k))) / (sqrt_2 * wL);
                    arg_x2_dis = (sqrt_pi * r_PD - 2 * ( - x_DE(k))) / (sqrt_2 * wL);
                    arg_y1_dis = (sqrt_pi * r_PD + 2 * ( - y_DE(k))) / (sqrt_2 * wL);
                    arg_y2_dis = (sqrt_pi * r_PD - 2 * (- y_DE(k))) / (sqrt_2 * wL);
                    term_x_dis = (erf(arg_x1_dis) + erf(arg_x2_dis));
                    term_y_dis = (erf(arg_y1_dis) + erf(arg_y2_dis));
                    Hii_MIMO_dis(i,i) = 0.25.*term_x_dis* term_y_dis;
                    H_MIMO_dis(i,i) = Hii_MIMO_dis(i,i);

            end
            end
        end

        var_dis(i) = (4*K_boltz*T/RL)*B*Fn + 2*q*(sum(R_PD*Hij_MIMO_dis(i,:)*Pt))*B + RIN*(sum((R_PD*Hij_MIMO_dis(i,:)*Pt).^2))*B;
        SNR_dis(i) = ((R_PD^2) * (Hii_MIMO_dis(i,i).^2) * P_elec) / (sum((R_PD^2).*(Hij_MIMO_dis(i,:).^2).*P_elec)+ var_dis(i));
        Rate_dis(i) = zeta*B*log2(1 + SNR_dis(i)/Gamma);  % Transmission rate for VCSEL_i

    end

    % SVD & Reassociation
    [~,S_dis,~] = svd(H_MIMO_dis);
    % Calculate the threshold based on the percentage of the maximum singular value
    threshold = 0.001 * max(S_dis(:));
    % Determine the number of significant singular values
    % kk = sum(S_dis > threshold);
    % kk = floor(sqrt(Nr(l))) - 1;
    kk = Nt;
    % Select the highest eigenvalues and corresponding singular vectors
    lambdai = diag(S_dis(:,1:sum(kk)));
    % Create a new diagonal matrix with the highest eigenvalues
    newS = diag(lambdai);
    % newS = S_dis;
    var_dis_svd = var_dis(1:sum(kk));
    SNR_svd_dis = ((R_PD^2)*P_elec.*(newS.^2)) ./ var_dis_svd;
    % SNR_svd_dis(i) = ((R_PD^2)*P_elec.*(S_dis(i,i).^2)) ./ var_dis(i);
    % SNR_svd_dis(i) = ((R_PD^2)*P_elec.*(H_MIMO_dis(i,i).^2)) ./ var_dis(i);
    % Rate_svd_dis(i) = zeta*B*log2(1 + SNR_svd_dis(i)/Gamma);  % Transmission rate for VCSEL_i
    Rate_svd_dis = zeta.*B.*log2(1 + SNR_svd_dis./Gamma);  % Transmission rate for VCSEL_i
    sum_Rate_dis(k) =  sum(Rate_dis(:))./1e12; % sum rate of MIMO OFDM (Tbps)
    sum_Rate_svd_dis(k) =  sum(Rate_svd_dis(:))./1e12; % sum rate of multiple identical SISO OFDM links (Tbps)

end
wosvd = sum_Rate_dis;
wsvd = sum_Rate_svd_dis;

% plot(r_DE.*1e3,sum_Rate_dis,'LineWidth',2)
% plot(r_DE.*1e3,sum_Rate_svd_dis,'LineWidth',2)


plot(r_DE.*1e3,wosvd,'b','LineWidth',2)
plot(r_DE.*1e3,wsvd,'r','LineWidth',2)
% plot(r_DE.*1e3,wsvd(2,:),'y','LineWidth',2)
% plot(r_DE.*1e3,wsvd(3,:),'m','LineWidth',2)