% c = 2.75; beta = 1.8e-9;
% 
% delta_I0 = 0.02; delta_q = 0.1;
% 
% gama_I = 0.012; gama_H = 0.02;
% 
% q = 10e-7; alpha = 0.1008;
% 
% theta = 1; lam = 1 / 14;
% 
% T = 40; t = 1; NN = T / t;
c=2.75;beta=1.74e-9;

delta_I0=0.1; delta_q=0.13;

gama_I=0.007; gama_H=0.02;

q=2e-7; alpha=1.2405e-4;

theta=0.75; lam=1/14;

T=60; t=0.01; NN=T/t;

% Initial values

% S = 5917 * 1e4; E = 4007; I = 524 * 1.5; Sq = 2776;
% 
% Eq = 400; H = I + Eq; R = 31; De = 24; sigma = 1 / 7;

 S = 60481283; E = 620; I = 33 * 2; Sq = 0;
 Eq = 2; H = I + Eq; R = 2; sigma = 1 / 7;
AA = [S E I Sq Eq H R];

for ii = 1:NN

% Increased isolation speed of patient due to new hospitals

if (ii * t) >= 14

delta_I = delta_I0 * 1.2;

else

delta_I = delta_I0;

end

% Modified SEIR Transmission dynamics model

dS = -(beta * c + c * q * (1 - beta)) * S * (I + theta * E) + lam * Sq;

dE = beta * c * (1 - q) * S * (I + theta * E) - sigma * E;

dI = sigma * E - (delta_I + alpha + gama_I) * I;

dSq = (1 - beta) * c * q * S * (I + theta * E) - lam * Sq;

dEq = beta * c * q * S * (I + theta * E) - delta_q * Eq;

dH = delta_I * I + delta_q * Eq - (alpha + gama_H) * H;

dR = gama_I * I + gama_H * H;

% Euler integration algorithm

S = S + dS * t;

E = E + dE * t;

I = I + dI * t;

Sq = Sq + dSq * t;

Eq = Eq + dEq * t;

H = H + dH * t;

R = R + dR * t;

AA = [AA; S E I Sq Eq H R];

end

% Theoretical estimation

Infected(:, 1) = round(AA(1:1 / t : size(AA, 1), 3));

Cured(:, 1) = round(AA(1:1 / t : size(AA, 1), 7));

plot(0:T, [Infected Cured])

% raw data

hold on

Data = importdata('意大利.csv');

data_Infected = Data.data(:, 1);

data_Cured = Data.data(:, 3);

loss1 = 0;
for i = 1:31
    loss1 = loss1 + (data_Infected(i,1) - Infected(i,1))^2;
end
loss1 = loss1/31

loss2 = 0;
for i = 1:31
    loss2 = loss2 + (data_Cured(i,1) - Cured(i,1))^2;
end
loss2 = loss2/31


plot([1:length(data_Infected)],[data_Infected data_Cured],' * ')
