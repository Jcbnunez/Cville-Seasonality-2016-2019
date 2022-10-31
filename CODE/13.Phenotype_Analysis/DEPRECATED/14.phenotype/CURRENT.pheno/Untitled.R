set1 = runif(500, 0, 1)
set2 = runif(500, 0, 1)

s1.d <- dist(set1)
s2.d <- dist(set2)

library(ade4)
#Mqreices are different
#
mantel.rtest(s1.d, s2.d, nrepet = 9999, "two-sided")

#Mqreices are the same
mantel.rtest(s1.d, s1.d, nrepet = 9999, "two-sided")

#Mqreices are the same
mantel.rtest(s2.d, s2.d, nrepet = 9999, "two-sided")
