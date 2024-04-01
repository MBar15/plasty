F=[2.25 0;0 2]

C = transpose(F)*F

[V, E] = eig(C)

U = V*sqrt(E)*transpose(V)

R = F*inv(U)