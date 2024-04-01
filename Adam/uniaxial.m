clc;
clear all;
close all;


syms c10 c20 lambda I1 sigma1 p sigma2 c01 c20 c11 c02 I2 sigma epsilonLn P lambda1 lambda2 lambda3

%lambda = 0.7

W = c10*(I1-3)+c20*(I1-3)^2

F = [lambda 0 0 ; 0 sqrt(1/lambda) 0; 0 0 sqrt(1/lambda)]

J = det(F) == 1
lambda2 = sqrt(1/lambda) == lambda3


B = F*transpose(F)

dWI1 = diff(W,I1)

sigmaT = [sigma 0 0; 0 0 0; 0 0 0] 

rce1 = sigmaT == -p*eye(3)+2*dWI1*B

p =solve(rce1(3,3),p)

I1 = trace(B)
c10 = 0.8;
c01 = 1;
c20 = 0.5;
c11 = 0.2;
c02 = 0.1;

p=subs(p);
rce1 = subs(rce1)

sigma1 = rhs(rce1(1,1))
sigma1 = subs(sigma1);

fplot(sigma1,[1 1.3],'color','green')
legend('sigma1')
xlabel("Protažení lambda")
ylabel("Napětí - PK")

sigma_biax_final = vpa(subs(sigma1,lambda,1.5),4)