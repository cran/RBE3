dgumbel2 <-
function(x,log=FALSE)
{lf<- x-exp(x)
if(log) lf
else exp(lf)
}
