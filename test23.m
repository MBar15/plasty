clc;
clear all;
close all;

% Vykreslete křivku 𝑃 − 𝜀𝑙𝑛 do stlačení vzorku 30% jednoosým tlakovým testem. Materiál je 
% nestlačitelný hyperelastický s modelem Mooney –Rivlin 5-ti parametry a konstantami: 
% c10=0.4MPa, c01=-1MPa, c20=0,2MPa, c11=0,2MPa, c02=0,1MPa. Jak se mění 20mm průměr 
% 15mm vysokého vzorku při tomto stlačení?

syms c10 c01 c11 c20 c02 I1 I2 lambda p


W = c10*(I1-3) + c01*(I2-3) +c11*(I1-3)*(I2-3) +c20*(I1-3)^2 + c02*(I1-3)^2
dWdI1 = diff(W,I1)
dWdI2 = diff(W,I2)

lambda1 = lambda
lambda2 = lambda1^-2
lambda3 = lambda1^-2

F = [
    lambda1,0,0;
    0,lambda2,0;
    0,0,lambda3]

B = F*transpose(F)

sigma = -p*eye(3) + 2*dWdI1*B - 2*dWdI2*inv(B) 

% stlačení sigmay a sigmaz je rovno 0 jen sigma x má hodnotu
rce1 = 0 == sigma(3,3)
p = solve(rce1,p)

% c10=0.4MPa, c01=-1MPa, c20=0,2MPa, c11=0,2MPa, c02=0,1MPa.
c10 = 0.4
c01 = -1
c11 = 0.2
c20 = 0.2
c02 = 0.1
I1 = lambda1^2 + lambda2^2 + lambda3^2
I2 = lambda1*lambda2 + lambda2*lambda3 + lambda3*lambda1

sigma = subs(sigma)
sigma = subs(sigma)

% prvni piolakirchoff
P = sigma(1,1)/lambda
lambda_val = linspace(0.7,1,50)
lambda_val_ln = log(lambda_val)
P_val = subs(P,lambda,lambda_val)
plot(lambda_val_ln,P_val)
xlim([lambda_val_ln(1) lambda_val_ln(end)])

