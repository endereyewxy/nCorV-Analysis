# 意大利新型冠状病毒发展预测

## 问题重述
新冠病毒已经在世界100多个国家蔓延，请选择一个国家，收集该国近一个多月的确诊病人数量，治愈数量和死亡人数，通过建立疾病传播模型，估计参数，预测其未来确诊病人数量，作出相应的图形，对结果进行分析。
## 问题分析
### 数据采集

### 传统模型
SEIR传播动力学模型
### 资料分析

## 模型假设

## 模型建立及求解

根据原理，可建立微分方程组：
$$
\left\{
\begin{split}
\frac{dS}{dt} &= -[c \beta + c q (1 - \beta)]S(I + \theta E) + \lambda S_q \\
\frac{dE}{dt} &= c \beta (1 - q) S (I + \theta E) - \sigma E \\
\frac{dI}{dt} &= \sigma E - (\delta_I + \alpha + \gamma_I) I \\
\frac{dS_q}{dt} &= c q (1 - \beta) S (I + \theta E) - \delta_q E_q \\
\frac{dE_q}{dt} &= c q \beta S (I + \theta E) - \delta_q E_q \\
\frac{dH}{dt} &= \delta_I I + \delta_q E_q - (\alpha + \gamma_H) H \\
\frac{dR}{dt} &= \gamma_I I + \gamma_H H \\
\frac{dD}{dt} &= \alpha(I + H)
\end{split}
\right.
$$

其中，$ \sigma $ 为潜伏着转化为感染者的速率， $ \alpha $ 为死亡率（同时作用于感染者和隔离治疗人员）， $ \delta_I $ 和 $ \delta_q $ 分别为感染者和潜伏者被隔离的速率， $ \gamma_I $ 和 $ \gamma_H $ 分别为感染者和隔离治疗人员康复的速率。

方便起见，该微分方程组不进行数值求解，而是以 $ t = 0.01 $ 天为间隔进行迭代，即：
$$
\left\{
\begin{split}
S_{i + 1} &= S_i + \frac{dS}{dt} \cdot t \\
E_{i + 1} &= E_i + \frac{dE}{dt} \cdot t \\
I_{i + 1} &= I_i + \frac{dI}{dt} \cdot t \\
S_{q(i + 1)} &= S_{qi} + \frac{dS_q}{dt} \cdot t \\
E_{q(i + 1)} &= E_{qi} + \frac{dE_q}{dt} \cdot t \\
H_{i + 1} &= H_i + \frac{dH}{dt} \cdot t \\
R_{i + 1} &= R_i + \frac{dR}{dt} \cdot t \\
D_{i + 1} &= D_i + \frac{dD}{dt} \cdot t \\
\end{split}
\right.
$$
其中 $ S_i $ 为第 $ i $ 次迭代的值， $ E_i, I_i, S_{qi}, E_{qi}, H_i, R_i, D_i $ 类似。

## 结果及分析

<img src="simulate.svg" />

## 参考文献

1. Biao Tang, Xia Wang, Qian Li, et al. Estimation of the transmission risk of 2019-nCov and its implication for public health interventions. 2020, Doi: 10.2139/ssrn.3525558

## 附录

```matlab
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
```