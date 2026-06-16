Nsymb = 2^15;
symb_rate = 100e6;
M_order = 2; % bits per symble
Vpi = 1;

P_in = 1e-3;
p_in_dBm = 10*log10(P_in);
P_Lo = 20e-3;
P_Lo_dBm = 10*log10(P_Lo);
E_in = sqrt(P_in);
E_Lo = sqrt(P_Lo);
linewidth = 1000;
T_symbol = 1/symb_rate;
ER_dB = 10*log10(1000);
IL_dB = 1.2;
Vbias = 1;
osnr_dB = 15;
R_resp = 0.7;
B_elec = symb_rate/2;
q = 1.6 * 10^(-19);

%generat random symbols

X1 = randi([0,1],2*Nsymb,1);
X2 = reshape(X1,2,[])';
X3 = X2 ;

I = 1 - 2*X3(:,2);
Q = 1 - 2*X3(:,1);
symbols_in = I + 1j*Q;
scatter(real(symbols_in), imag(symbols_in), 'filled')
grid on
axis equal

xlabel('I')
ylabel('Q')
title('QPSK Constellation')

% voltage generation using the symbols 
V_I = I*Vpi/2;
V_Q = Q*Vpi/2;

% RIN add
RIN_dB = -155;   % dB/Hz
RIN_linear = 10^(RIN_dB/10);
B_in = symb_rate/2;
sigma_rin = sqrt(RIN_linear*B_in);

power_noise_in = sigma_rin*randn(Nsymb,1);
power_noise_LO = sigma_rin*randn(Nsymb,1);

E_in_with_noise = sqrt(P_in.*(1 + power_noise_in));

gamma = 0.5*sqrt(1 - ((10^(ER_dB/10)-1)/(10^(ER_dB/10)+1))^2);

phi_I = pi*(Vbias/Vpi + V_I/Vpi);
phi_Q = pi*(Vbias/Vpi + V_Q/Vpi);

const_E_in = E_in_with_noise./(2*10^(IL_dB/20));

E_out_I = const_E_in .* ...
    (sqrt(1+2*gamma).*exp(1j*phi_I/2) + ...
     sqrt(1-2*gamma).*exp(-1j*phi_I/2));

E_out_Q = 1j*const_E_in .* ...
    (sqrt(1+2*gamma).*exp(1j*phi_Q/2) + ...
     sqrt(1-2*gamma).*exp(-1j*phi_Q/2));

E_total = E_out_I + E_out_Q;

p_signal = mean(abs(E_total).^2);
OSNR_lin = 10^(osnr_dB/10);
P_noise = p_signal/OSNR_lin;

noise_ase = sqrt(P_noise/2) .* ...
    (randn(Nsymb,1) + 1j*randn(Nsymb,1));

E_after = E_total + noise_ase;

sigma2_S = 2*pi*linewidth*T_symbol;
sigma2_LO = 2*pi*linewidth*T_symbol;

phi_S = cumsum(sqrt(sigma2_S).*randn(Nsymb,1));
phi_LO = cumsum(sqrt(sigma2_LO).*randn(Nsymb,1));

E_s = E_after .* exp(1j*phi_S);
E_LO = sqrt(P_Lo.*(1 + power_noise_LO)) .* exp(1j*phi_LO);

E1 = (E_s + E_LO)/sqrt(2);
E2 = (E_s - E_LO)/sqrt(2);
E3 = (E_s + 1j*E_LO)/sqrt(2);
E4 = (E_s - 1j*E_LO)/sqrt(2);

I1_mean = R_resp*abs(E1).^2;
I2_mean = R_resp*abs(E2).^2;
I3_mean = R_resp*abs(E3).^2;
I4_mean = R_resp*abs(E4).^2;

sigma1 = sqrt(2*q*B_elec.*I1_mean);
sigma2 = sqrt(2*q*B_elec.*I2_mean);
sigma3 = sqrt(2*q*B_elec.*I3_mean);
sigma4 = sqrt(2*q*B_elec.*I4_mean);

i1 = I1_mean + sigma1.*randn(Nsymb,1);
i2 = I2_mean + sigma2.*randn(Nsymb,1);
i3 = I3_mean + sigma3.*randn(Nsymb,1);
i4 = I4_mean + sigma4.*randn(Nsymb,1);

i_x = i1 - i2;
i_p = i3 - i4;

x_B = i_x ./ (2*R_resp.*abs(E_LO));
p_B = i_p ./ (2*R_resp.*abs(E_LO));

figure
scatter(real(E_total),imag(E_total),'filled')
grid on
axis equal

figure
scatter(x_B,p_B,'filled')
grid on
axis equal

% First-Order Coherence Function
max_lag = 20000;
tau = max_lag*T_symbol;
lags = (0:max_lag).';
tau = lags*T_symbol;
theory_gamma = exp(-pi*linewidth*tau);
tau_sim = lags*T_symbol;
gamma_sim = zeros(max_lag+1,1);
E_field = exp(1j*phi_S);

for k = 0:max_lag
    gamma_sim(k+1) = abs(mean(conj(E_field(1:end-k)).*E_field(1+k:end)));
end


figure
plot(tau_sim*1e6,gamma_sim,'LineWidth',1.5)
hold on
plot(tau*1e6,theory_gamma,'--','LineWidth',1.5)
grid on

xlabel('\tau (\mus)')
ylabel('|g^{(1)}(\tau)|')
title('First-Order Coherence Function')
legend('Simulation','Theory')