clc;
clear all;
close all;


syms c10 c20 lambda I1 sigma1 p sigma2 c01 c20 c11 c02 I2 sigma epsilonLn P

%lambda = 0.7

W = c10*(I1-3)+c01*(I2-3)+c11*(I1-3)*(I2-3)+c20*(I1-3)^2+c02*(I2-3)^2

F = [lambda^(-1/2) 0 0 ; 0 lambda^(-1/2) 0; 0 0 lambda]

J = det(F)

B = F*transpose(F)

dWI1 = diff(W,I1)
dWI2 = diff(W,I2)

sigmaT = [sigma 0 0; 0 sigma 0; 0 0 0] 

rce1 = sigmaT == -p*eye(3)+2*dWI1*B-2*dWI2*inv(B)

p =solve(rce1(3,3),p)

I1 = trace(B)
I2 = 1/2*(trace(B)^2 - trace(B^2))
c10 = 0.4;
c01 = 1;
c20 = 0.2;
c11 = 0.2;
c02 = 0.1;
d= 20;
p=subs(p);
rce1 = subs(rce1)

sigma1 = rhs(rce1(1,1))

% Převod: 
PK11=sigma1/lambda

fplot(PK11,[log(0.7) log(1)],'color','green')
legend('FPK11')
xlabel("Protažení lambda")
ylabel("Napětí - PK")

sigma_biax_final = vpa(subs(sigma1,lambda,0.7),4)
PK_biax_final = vpa(subs(PK11,lambda,0.7),4)

%Změna průměru

% F = PK11*3.1415*d^2/4
% 
% d1 = sqrt(4*F/(3.1415*sigma1))
% 
% fplot(d1,[0.7 1])






















