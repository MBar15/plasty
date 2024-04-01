clear all
close all
clc
%zadaní. Urcete FPK napětí v obou směrech v materiálu pro biaxiální zkoušku s poměrem
%protažení 2:1 pro protaženi lambdaX=2.4 lambdaY=1.7. Materiál je definován
%Yeoh modelem 3. řádu s konstantami c10=0.6MPa, c20=-0.1MPa, c30=0.02MPa a
%je objemově nestlačitelný

%nejprve si nadefinujeme potřebné parametry
syms c10 c20 c30 I1 lambda11 lambda22 lambda33 p lambda

%známe funkci energie napjatosti
W=c10*(I1-3)+c20*(I1-3)^2+c30*(I1-3)^3
%potřebujeme její derivaci podle I1
dW=diff(W,I1)
%Tenzor deformačního gradientu. Jedná se o namáhání v normálovém směru,
%takže mimodiagonální složky tam nejsou a diagonáljí jsou přímo rovny
%hlavním protažením
F=[lambda11 0 0;0 lambda22 0; 0 0 lambda33]
%levý Cauchy-Greenův tenzor deformace
B=F*transpose(F)
%ted potřebujeme nadefinovat, čemu se vlastně I1 rovná
I1=lambda11^2+lambda22^2+lambda33^2
%také je třeba použít podmínku obj. nestlačitelnosti, tedy, že jakobián se
%rovná 1
J=lambda11*lambda22*lambda33==1
%a také nadefinovat jaké jsou poměry lambd při dvouosé zkoušce kdy jeden
%směr tahám 2x více než druhý. Nelze prostě lambdu11 dělit dvěma, protože
%nezatížený stav je lambda=1.
lambda22=1+(lambda11-1)/2
% a protože už známe poměr protažení, tak lze vše vyjadřovat jako jednu lambdu bez indexů
lambda11=lambda
%přeznačení dokončíme dosazením do lambda22
lambda22=subs(lambda22)
% nyní lze z výše uvedené rovnice pro jakobián vyjádřit třetí lambdu 
lambda33=solve(J,lambda33)
%dosadíme výsledek, do I1
I1=subs(I1)
%i do derivace W
dW=subs(dW)
%i do B
B=subs(B)
%jednotkový tenzor je definován takto:
I=[1 0 0;0 1 0;0 0 1]

%známý obecný vztah pro výpočet Cauchyho napětí z energie napjatosti
sigma=2*dW*B-p*I
%protože se jedná o nestlačitelný materiál, tak nelze napětí dopočítat,
%protože ve výše uvedené rovnici neznáme hodnotu hydrostatického tlaku p.
%Použijeme tedy okrajovou podmínku. Při dvouosé zkoušce je napětí ve třetím
%směru nulové. Tedy:
rce1=sigma(3,3)==0
%Z této rovnice lze dopočítat chybějící hydrostatický tlak
p=solve(rce1,p)
%a ten nyní dosadíme ro dovnice pro napětí
sigma=subs(sigma)

%nyní už můžeme zadat materiálové konstanty materiálu
c10=0.6
c20=-0.1
c30=0.02


%vše lze dosadit do rovnice pro Cauchyho napětí
sigma=subs(sigma)
%v zadání je však, že chceme výsledky v inženýrských napětích FPK. Je tedy
%nutné je převést na FPK. První hlavní napětí dělíme první hlavní lambdou a
%dostaneme první hlavní napětí FPK
FPK11=sigma(1,1)/lambda11
%analogicky pro FPK22
FPK22=sigma(2,2)/lambda22
%nyní už můžeme obě napětí vykreslit v požadovaných rozsazích
fplot(FPK11,[1 2.4],'color','green')
hold on
fplot(FPK22,[1 1.7],'color','red')
legend('FPK11','FPK22')

%na závěr dopočítáme konkrétní číselné hodnoty v MPa pro lambda11=2.4 a lambda22=1.7
%příkaz vpa slouží k numerickému vyčíslení
FPK11_end=vpa(subs(FPK11,lambda11,2.4),3)
FPK22_end=vpa(subs(FPK22,lambda11,1.7),3)

%Pro zadaný materiál, tyzp zkoušky a rozsah protažení jsou hlavní inženýrské
%napětí 6.7MPa  a 1.05MPa.