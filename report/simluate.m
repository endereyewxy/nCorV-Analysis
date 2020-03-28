hold on

% Model parameters

c        = 2.75;
beta     = 1.74e-9;
delta_I0 = 0.1;
delta_q  = 0.13;
gama_I   = 0.001;
gama_H   = 0.03;
q        = 2e-7;
alpha    = 0.009;
theta    = 0.75;
lam      = 1/14;
sigma    = 1 / 7;

% Simulation parameters

T  = 60;
t  = 0.01;
NN = T/t;

% Initial values

S  = 60481283;
E  = 620;
I  = 33 * 2;
Sq = 0;
Eq = 2;
H  = I + Eq;
R = 2;
D = 0;

G = zeros(NN, 8);
G(1, :) = [S E I Sq Eq H R D];

for ii = 1 : NN
    if (ii * t) >= 23
        delta_I = delta_I0 * 1.7;
    else
        delta_I = delta_I0;
    end
    %delta_I = delta_I0;
    
    dS  = -(beta * c + c * q * (1 - beta)) * S * (I + theta * E) + lam * Sq;
    dE  = beta * c * (1 - q) * S * (I + theta * E) - sigma * E;
    dI  = sigma * E - (delta_I + alpha + gama_I) * I;
    dSq = (1 - beta) * c * q * S * (I + theta * E) - lam * Sq;
    dEq = beta * c * q * S * (I + theta * E) - delta_q * Eq;
    dH  = delta_I * I + delta_q * Eq - (alpha + gama_H) * H;
    dR  = gama_I * I + gama_H * H;
    dD  = alpha * I + alpha * H;

    S  = S  + dS  * t;
    E  = E  + dE  * t;
    I  = I  + dI  * t;
    Sq = Sq + dSq * t;
    Eq = Eq + dEq * t;
    H  = H  + dH  * t;
    R  = R  + dR  * t;
    D  = D  + dD  * t;
    
    G(ii + 1, :) = [S E I Sq Eq H R D];
end

% Simulation result

yI(:, 1) = round(G(1 : 1/t : size(G, 1), 3));
yR(:, 1) = round(G(1 : 1/t : size(G, 1), 7));
yD(:, 1) = round(G(1 : 1/t : size(G, 1), 8));

plot(0 : T, [yI, yR, yD])

% Raw data

raw = importdata('意大利.csv');

rI = raw.data(:, 1);
rR = raw.data(:, 3);
rD = raw.data(:, 4);

plot(1 : length(rI), [rI, rR, rD],'x');

legend 预计感染 预计治愈 预计死亡 实际感染 实际治愈 实际死亡
title 意大利疫情预测
xlabel 天
ylabel 人
