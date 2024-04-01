clear all
close all
clc

syms lambda2 N11 N12 N21 N22 lambda3
%define deformation gradient tensor
F=[2.25 0;0 2]
%calculate right Cauchy green deformation tensor
C=transpose(F)*F
%identity tensor
I=[1 0;0 1]
%calculation of eigen stretches Alternative function is eig(C)
q=C-lambda2*I
q1=det(q)==0
lambda2=solve(q1,lambda2)

% sum of out of diagonal values in F. it is necessary to avoid errors in
% cases where there are no out of diagonal values in F. Because then it is
% not clear wheter orientation of eigen vector should be (1,0) or (-1,0).
% Further it is necessary to properly rank eigen values in such case

sumFij=F(1,2)+F(2,1)
if sumFij==0
    if F(1,1)>F(2,2)
       lambda11=vpa(sqrt(max(lambda2)),5)
       lambda22=vpa(sqrt(min(lambda2)),5)
    else
       lambda22=vpa(sqrt(max(lambda2)),5)
       lambda11=vpa(sqrt(min(lambda2)),5) 
    end
else
       lambda11=vpa(sqrt(max(lambda2)),5)
       lambda22=vpa(sqrt(min(lambda2)),5)    
end
% definition of eigen vectors
N1=[N11;N12]
N2=[N21;N22]
%calculation of first eigen vector. matrix equation
rce0=(C-lambda11^2*I)*N1==[0;0]
%eigen vector has unit length
rce1=N12==sqrt(1-(N11^2))
%take second equation from matrix equation. It does not matter whic eq. is taken
%if there are non diagonal terms in F
rce2=subs(rce0(2,1))
%express second coordinate of first egiven vector
N12=solve(rce2,N12)
%substitute to eq. describing eigen vector unit length
rce1=subs(rce1)
%solve first coordinate of eigen vector from eq. describing eigen vector
%unit length
N11=solve(rce1,N11)
%condition to account for cases with zeros in out of diagonal terms of F
%then positive root should be taken. Otherwise eq. rce1 has only one root
if sumFij==0
N11 = N11(N11>=0)
end
%substitute obtained values
N12=subs(N12)
% final first eigen vector
N1=subs(N1)

%the same operations are done for second eigenvector
rce3=(C-lambda22^2*I)*N2==[0;0]

rce4=N22==sqrt(1-(N21^2))
rce5=subs(rce3(2,1))
N21=solve(rce5,N21)

rce4=subs(rce4)

N22=solve(rce4,N22)
if sumFij==0
N22 = N22(N22>=0)
end
N21=subs(N21)
N2=subs(N2)
%stretch tensor caluclation
U=lambda11*N1*transpose(N1)+lambda22*N2*transpose(N2)
%rotation tensor calculation
R=F*inv(U)

%Extension 2D deformation in 3D for incompressible material
U3=[U(1,1) U(1,2) 0;U(2,1) U(2,2) 0;0 0 1/(U(1,1)*U(2,2))]
%calculation of principal stretches
lambdaP=eig(U3)