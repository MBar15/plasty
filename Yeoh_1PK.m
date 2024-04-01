clear all
close all
clc
%zadaní. Urcete FPK napětí v obou směrech v materiálu pro biaxiální zkoušku s poměrem
%protažení 2:1 pro protaženi lambdaX=2.4 lambdaY=1.7. Materiál je definován
%Yeoh modelem 3. řádu s konstantami c10=0.6MPa, c20=-0.1MPa, c30=0.02MPa a
%je objemově nestlačitelný

syms lambda c10 c20 c30 I1 p

W = c10*(I1-3) + c20*(I1-3)^2 + c30*(I1-3)^3
dW = diff(W,I1)

lambda1= lambda
lambda2 = (lambda1-1)/2 + 1
lambda3 = 1/(lambda1*lambda2)
F = [lambda1 0 0;
     0 lambda2 0;
     0 0 lambda3]

B = F*transpose(F)
I1 = lambda1^2 + lambda2^2 + lambda3^2 

c10 = 0.6;
c20 = -0.1;
c30 = 0.02;

sigma = -p*eye(3) + 2*dW*B
sigma = subs(sigma)

rce1 = 0 == sigma(3,3)
p = solve(rce1,p)
sigma = subs(sigma)

FPK11 = sigma(1,1)/lambda1
FPK22 = sigma(2,2)/lambda2

fplot(FPK11,[1 2.4],'color','green')
hold on
fplot(FPK22,[1 2.4],'color','red')
legend('FPK11','FPK22')


% lamda je jen jedna takže napětí jsou funkce lamda1 = lamda a tudíž musím
% do obou rovnic dosadit 2.4 protože napětí 2 je funkcí lamdy1 jako
% polovina protažení, si myslím že přehlédl chybu
FPK11_end=vpa(subs(FPK11,lambda,2.4),3)
FPK22_end=vpa(subs(FPK22,lambda,2.4),3)