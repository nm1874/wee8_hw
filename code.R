#1. 

#part a

#let v = x'
#w <- c(x,v)
#w_ <- A%*%w we must find A 
#w_ <- c(v,x'')
#x'' = -3v -2x or -2x -3v
#x'' = c(-2, -3) %*% c(x, v)

#v = 0x + 1v
#v = c(0,1)%*% c(x, v)

A <- matrix(c(0, -2, 1, -3), 2); A

#part b: find exp(At)

P <- eigen(A)$vectors; P
PInv <- solve(P);
eigen(A)$values


Aexp <- function(t) P%*%diag(c(exp(-2*t),exp(-1*t)))%*%PInv

#part c: 
par(mar=c(4,1,1,1)) #leave room at bottom for the equation
plot(NULL, xlim = c(0,4),ylim = c(0,4), xlab = "", ylab = "", axes = FALSE, 
     asp = 1, pch = 20)
axis(1, pos = 0); axis(2, pos = 0)
mtext(paste("x.dot  =", A[1,1],"x + ", A[1,2], "y"),1,1)
mtext(paste("y.dot  =", A[2,1],"x + ", A[2,2], "y"),1,2)

#Choose an initial value for the vector v
v0 <- c(1,0)

#Choose a sequence of times at which to plot the solution
times <- seq(from = 0, to =4, by = 0.1)
for (t in times)
  points(t, (Aexp(t)%*%v0)[1], pch = 20 )
v0 <- c(0,1)
for (t in times)
  points(t, (Aexp(t)%*%v0)[1], pch = 20, col = "green" )
#If v0 is an eigenvector the solution moves out or in along a line
v0 <- c(1,-3)
for (t in times)
  points(t, (Aexp(t)%*%v0)[1], pch = 20, col = "red" )

#2. 
S<-matrix(c(1,2,2,1),2); S
Sexp <- function(t) exp(t)*matrix(c(cosh(2*t), sinh(2*t), sinh(2*t), cosh(2*t)),2)

plot(NULL, xlim = c(-4,4),ylim = c(-4,4), xlab = "", ylab = "", axes = FALSE, asp = 1, pch = 20)
axis(1, pos = 0); axis(2, pos = 0)
mtext(paste("x.dot  =", S[1,1],"x + ", S[1,2], "y"),1,1)
mtext(paste("y.dot  =", S[2,1],"x + ", S[2,2], "y"),1,2)

#Choose an initial value for the vector v
v0 <- c(1,0)
#Choose a sequence of times at which to plot the solution
times <- seq(from = 0, to =2, by = 0.1)
for (t in times)
  points((Sexp(t)%*%v0)[1], (Sexp(t)%*%v0)[2], pch = 20 )
v0 <- c(0,1)
for (t in times)
  points((Sexp(t)%*%v0)[1], (Sexp(t)%*%v0)[2], pch = 20, col = "green" )
v0 <- c(-1,0)
for (t in times)
  points((Sexp(t)%*%v0)[1], (Sexp(t)%*%v0)[2], pch = 20, col = "blue" )

#The matrix S defines a vector field that specifies
#the velocity vector as a function of position.
plotvec = function(p) {
  velocity <- S%*%p
  arrows(p[1],p[2],p[1]+ 0.1*velocity[1], p[2]+ 0.1*velocity[2], length = 0.1)
}

#We can draw a map of the vector field.
#Every vector is tangent to one of the solution curves.
for (x in -5:5)
  for (y in -3:3)
    plotvec(c(x,y))


#3. a)
B<-matrix(c(-1,-1,9,5),2); B
eigen(B)$values
L<- 2
v<-eigen(B)$vectors[,1]; v
N <- B - diag(c(L,L)); N; N%*%N

B<- diag(c(L,L)) + N; B; B%*%B

#part b
B<- diag(c(L,L)) + N;B
L
#exp(At)=exp(2It)%*%exp(Nt)
#matrices that commute 

#exp(Nt) = I + Nt + N^2t^2/2! + ... 
#since N is nilpotent, every exponent greater than for equal to 2 = the zero matrix 
#so exp(Nt) = I + matrix(c(-3t, -t, 9t, 3t),2)
Bexp <- function(t) diag(c(exp(L*t),exp(L*t)))%*%(diag(c(1,1))+ t*N)

plot(NULL, xlim = c(-4,4),ylim = c(-4,4), xlab = "", ylab = "", axes = FALSE, asp = 1, pch = 20)
axis(1, pos = 0); axis(2, pos = 0)
mtext(paste("x.dot  =", B[1,1],"x + ", B[1,2], "y"),1,1)
mtext(paste("y.dot  =", B[2,1],"x + ", B[2,2], "y"),1,2)

#Choose an initial value for the vector v
v0 <- c(1,0)
#Choose a sequence of times at which to plot the solution
times <- seq(from = 0, to =2, by = 0.1)
for (t in times)
  points((Bexp(t)%*%v0)[1], (Bexp(t)%*%v0)[2], pch = 20 )
v0 <- c(0,1)
for (t in times)
  points((Bexp(t)%*%v0)[1], (Bexp(t)%*%v0)[2], pch = 20, col = "green" )

#The matrix B defines a vector field that specifies
#the velocity vector as a function of position.
plotvec = function(p) {
  velocity <- B%*%p
  arrows(p[1],p[2],p[1]+ 0.1*velocity[1], p[2]+ 0.1*velocity[2], length = 0.1)
}

#We can draw a map of the vector field.
#Every vector is tangent to one of the solution curves.
for (x in -5:5)
  for (y in -3:3)
    plotvec(c(x,y))
