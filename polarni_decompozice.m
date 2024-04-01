% polální dekompozice
% https://www.continuummechanics.org/polardecomposition.html

F = [
    1, 0.495, 0.5;
    -0.333, 1,- 0.247;
    0.959, 0, 1.5
    ]

C = transpose(F)*F

[V, E] = eig(C)

U = V*sqrt(E)*transpose(V)

R = F*inv(U)