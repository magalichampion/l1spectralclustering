library(Matrix)

N1=2
N2=5
k=sample(N1:N2,1)
T1=3
T2=10
taille=sort(sample(T1:T2,k,replace=FALSE))

A=matrix(1,taille[1],taille[1])

for (i in 2:length(taille)){
A=bdiag(A,matrix(1,taille[i],taille[i]));
}

n=dim(A)[1]
A=A-diag(dim(A)[1])

p=0.5
B=matrix(rbinom(n*n,1,p),n,n)
B[lower.tri(B)] <- 0
B=B+t(B)
B=B-diag(diag(B))
A_perturbed=(A+B)-2 * floor((A+B)/2)
