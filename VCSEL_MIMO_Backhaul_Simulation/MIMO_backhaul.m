clc
clear
close all

% Tbps VCSEL paper simulation (direct link simulation) 
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

for k = 1:length(K)
    for l = 1:length(w0)
          H_MIMO = []; 
            for i = 1:K(k)      
                    K_sq = sqrt(K(k));
                    mi = floor(i/K_sq); % row in H matrix
                    ni = i - (floor(i/K_sq) - 1)*K_sq; % column in H matrix
                    xi = (- (K_sq - 1)/2 + ni - 1)*d_PD;
                    yi = ( (K_sq - 1)/2 + mi + 1)*d_PD;      
                            for j = 1:K(k)
                                mj = floor(j/K_sq); % row in H matrix
                                nj = j - (floor(j/K_sq) - 1)*K_sq; % column in H matrix
                                xj = (- (K_sq - 1)/2 + nj - 1)*d_PD;
                                yj = ( (K_sq - 1)/2 + mj + 1)*d_PD;
                                            if i ~= j
                                                integrand_ij = @(x, y) ((2 / (pi.*wL_sq(l)))) .* exp(-2 * ((x - xi + xj).^2 + (y - yi + yj).^2) ./ (wL_sq(l)));
                                                % Hij:
                                                Hij_MIMO(i,j) = integral2(integrand_ij, y_min, y_max, x_min, x_max);
                                                H_MIMO(i,j) = Hij_MIMO(i,j);
                                            else if i == j
                                                    integrand_ii = @(x, y) ((2 / (pi.*wL_sq(l)))) .* exp(-2 * ((x).^2 + (y).^2) ./ (wL_sq(l)));
                                                    % Hii 
                                                    Hii_MIMO(i,j) = integral2(integrand_ii, y_min, y_max, x_min, x_max);
                                                    H_MIMO(i,j) = Hii_MIMO(i,j);
                                              end
                                            end
                                        end
           var(i) = (4*K_boltz*T/RL)*B*Fn + 2*q*(sum(R_PD*Hij_MIMO(i,:)*Pt))*B + RIN*(sum((R_PD*Hij_MIMO(i,:)*Pt).^2))*B;
           SNR(i) = ((R_PD^2) * (Hii_MIMO(i,i).^2) * P_elec) / (sum((R_PD^2).*(Hij_MIMO(i,:).^2).*P_elec)+ var(i));
            Rate(i) = zeta*B*log2(1 + SNR(i)/Gamma);  % Transmission rate for VCSEL_i

            end   
           
            [~,S,~] = svd(H_MIMO);
           threshold = 0.01 * max(S(:));
           % Determine the number of significant singular values
           kk = sum(S > threshold);
           % Select the highest eigenvalues and corresponding singular vectors
           lambdai = diag(S(:,1:sum(kk)));
           % Create a new diagonal matrix with the highest eigenvalues
           newS = diag(lambdai);
           var_new = var(1:sum(kk));
           SNR_svd = ((R_PD^2)*P_elec.*(newS.^2)) ./ var(i);
            % SNR_svd(i) = ((R_PD^2)*P_elec.*(H_MIMO(i,i).^2)) ./ var(i);
            Rate_svd = zeta*B*log2(1 + SNR_svd/Gamma);  % Transmission rate for VCSEL_i

            sum_Rate(l) =  sum(Rate(:))/1e12; % sum rate of MIMO OFDM (Tbps)   
            sum_Rate_svd(l) =  sum(Rate_svd(:))/1e12; % sum rate of multiple identical SISO OFDM links (Tbps)  

    end

    plot(w0*1e6,sum_Rate,'b','LineWidth',2)
    plot(w0*1e6,sum_Rate_svd,'-.r','LineWidth',2)
    % plot(sqrt(wL_sq)*1e3,sum_Rate,'b','LineWidth',2)
    % plot(sqrt(wL_sq)*1e3, sum_Rate_svd,'-.r','LineWidth',2)
end

% xlabel('\omega(L) [mm]');

hold off 