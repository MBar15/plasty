%zadaní. Urcete FPK napětí v obou směrech v materiálu pro biaxiální zkoušku s poměrem
%protažení 2:1 pro protaženi lambdaX=2.4 lambdaY=1.7. Materiál je definován
%Yeoh modelem 3. řádu s konstantami c10=0.6MPa, c20=-0.1MPa, c30=0.02MPa a
%je objemově nestlačitelný


clc;
clear all;
close all;


syms c10 c20 lambda I1 sigma1 p sigma2 c01 c20 c11 c02 I2 sigma epsilonLn P c30 lambda1 lambda2 lambda3

W = c10*(I1-3)+c20*(I1-3)^2+c30*(I1-3)^3
dWI1 = diff(W,I1)

F = [lambda1 0 0 ; 0 lambda2 0; 0 0 lambda3]

J = det(F) == 1
lambda2 = 1/2 + lambda1/2
lambda1 = lambda 
lambda2 = subs(lambda2)

lambda3 = solve(J,lambda3)

B = F*transpose(F)

sigmaT = [sigma 0 0; 0 sigma 0; 0 0 0] 

rce1 = sigmaT == -p*eye(3)+2*dWI1*B

p =solve(rce1(3,3),p)

I1 = trace(B)
c10 = 0.6;
c01 = 1;
c20 = -0.1;
c11 = 0.2;
c02 = 0.1;
c30 = 0.02;
p=subs(p);
rce1 = subs(rce1)

sigma1 = rhs(rce1(1,1))
sigma2 = rhs(rce1(2,2))

sigma1 = subs(sigma1)
sigma2 = subs(sigma2)

% Převod: 
PK11=sigma1/lambda1
PK11=subs(PK11)
PK22=sigma2/lambda2
PK22=subs(PK22)

fplot(PK11,[1 2.4],'color','green')
hold on
fplot(PK22,[1 1.7],'color','black')
legend('FPK11','FPK22')
xlabel("Protažení lambda")
ylabel("Napětí - PK")

PK_biax_final_1 = vpa(subs(PK11,lambda,2.4),4)
PK_biax_final_2 = vpa(subs(PK22,lambda,1.7),4)























