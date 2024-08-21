NonIter1GFit<-function(x,y){

# non-iterative estimation, ised for parameter initialisation

x <- x[!is.nan(y)]; y <- y[!is.nan(y)];
x <- x[is.finite(y)]; y <- y[is.finite(y)];
x <- x^2;
y <- y-1;
n <- length(x);

# minimize absolute deviation
S = rep(0,n);
S[1]= 0;
S[2:n] = 1/2*(y[2:n]+y[1:n-1])*(x[2:n]-x[1:n-1]);
S = cumsum(S);

param <- solve ( matrix( c(sum((x-x[1])^2), sum((x-x[1])*S), sum((x-x[1])*S), sum(S^2) ), nrow=2, ncol=2, byrow = TRUE), 
                  matrix( c(sum((y-y[1])*(x-x[1])), sum(S*(y-y[1]))), nrow=2,  ncol=1, byrow = TRUE))
a1 = param[1];
c1 = param[2];


# b1, c1 - the weights of exponentials should be positive!

theta = exp(c1*x);
param2 <- solve ( matrix( c(length(x), sum(theta), sum(theta), sum(theta^2) ), # the data elements 
                        nrow=2,              # number of rows 
                        ncol=2,              # number of columns 
                        byrow = TRUE), 
                 matrix( c(sum(y), sum(theta*y)), nrow=2,  ncol=1, byrow = TRUE));
a1 = param2[1];
b1 = param2[2];


s = sqrt(-1/(4*c1));
kappa = 1/(b1*((4*pi)*s^2))^(3/2)

return( list("s" = s, "kappa" = kappa))
}