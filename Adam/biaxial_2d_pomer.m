clc;
clear all;
close all;


syms c10 c20 lambda I1 sigma1 p sigma2 c01 c20 c11 c02 I2 sigma epsilonLn P c30 lambda1 lambda2 lambda3

W = c10*(I1-3)+c20*(I1-3)^2
dWI1 = diff(W,I1)

F = [1.5*lambda 0 0 ; 0 1.25*lambda 0; 0 0 1/(1.5*1.25*lambda^2)]
F1 = subs(F,lambda,1)

J = det(F) == 1




B = F*transpose(F)

sigmaT = [sigma 0 0; 0 sigma 0; 0 0 0] 

rce1 = sigmaT == -p*eye(3)+2*dWI1*B

p =solve(rce1(3,3),p)

I1 = trace(B)
c10 = 0.8; 
c20 = 0.5;
p=subs(p);

rce1 = subs(rce1)

sigma1 = rhs(rce1(1,1))
sigma2 = rhs(rce1(2,2))

sigma1 = subs(sigma1)
sigma2 = subs(sigma2)

sigma1_biax_final = vpa(subs(sigma1,lambda,1),4)
sigma2_biax_final = vpa(subs(sigma2,lambda,1),4)

fplot(sigma1,[1 1.5])
hold on 
fplot(sigma2,[1 1.25],'g')

% Převod: 
PK11=sigma1/lambda
PK11=subs(PK11)
PK22=sigma2/lambda
PK22=subs(PK22)

fplot(PK11,[1 1.5],'color','green')
hold on
fplot(PK22,[1 1.25],'color','black')
legend('FPK11','FPK22')
xlabel("Protažení lambda")
ylabel("Napětí - PK")

PK_biax_final_1 = vpa(subs(PK11,lambda,1.5),4)
PK_biax_final_2 = vpa(subs(PK22,lambda,1.25),4)
