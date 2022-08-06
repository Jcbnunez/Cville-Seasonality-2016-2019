#by ABC
n.max = 16046
n.min = 1086

y = c(rep(n.max, 13), rep(n.min, 2))
h.m.abc = 1/mean(1/y)

#by SSE
n.max = 10000
n.min = 3000

y = c(rep(n.max, 13), rep(n.min, 2))
h.m.sse = 1/mean(1/y)