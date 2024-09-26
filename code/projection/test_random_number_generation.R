

a = array(NA, dim = c(2,3,4))


for (i in 1:2) {
  # set.seed(1234)
  for (j in 1:3) {
    set.seed(1234)
    for (k in 1:4) {
      # set.seed(1234)
      
      a[i,j,k] = sample(1:10, 1)
      
    }
    
  }
}

a[1,1,]
a[1,2,]