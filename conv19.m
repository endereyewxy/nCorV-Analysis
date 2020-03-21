% SEIR Transmission dynamics model of 2019 nCoV coronavirus with considering the weak infectious ability and changes in latency duration % %%%%%

%Shi Pengpeng, Cao Shengli, Feng Peihua, Completed on Jan 5th, first released on Jan 7th % %

%This latest modified version was released on Feb. 10th, 2020 % %

%Model parameters

c = 2.6; beta = 1.74e-9;

delta_I0 = 0.13; delta_q = 0.13;

gama_I = 0.0035; gama_H = 0.007;

q = 10e-7; alpha = 0.0001;

theta = 1; lam = 1 / 14;

T = 40; t = 0.01; NN = T / t;

% Initial values

S = 5917 * 1e4; E = 4007; I = 524 * 1.5; Sq = 2776;

Eq = 400; H = I + Eq; R = 31; De = 24; sigma = 1 / 7;

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

data_Infected = [524 658 958 1303 2567 3349 4334 5486 6738 8565 10532 12712 15679 18445 20677 23139 24881 26965 28532 29659]';

data_Cured = [31 32 42 44 47 80 90 116 166 215 295 396 520 671 817 1115 1439 1795 2222 2639]';

plot([1:length(data_Infected)]'-1,[data_Infected data_Cured],' * ')