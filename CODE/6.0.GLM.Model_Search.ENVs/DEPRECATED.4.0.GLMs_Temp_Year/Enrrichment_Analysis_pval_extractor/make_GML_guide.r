### Make guide file
### 

pops <- c(
  "VA_ch",
  "UA_Ode",
  "TR_Yes",
  "PA_li",
  "FI_Aka",
  "DE_Mun",
  "DE_Bro"
)

p_vals <-c(
  seq(from = 0.1, to = 1, by = 0.1)[-10],
  seq(from = 0.01, to = 0.1, by = 0.01)[-10],
  seq(from = 0.001, to = 0.01, by = 0.001)[-10]
  ) %>% .[order(.)]

expand.grid(pops,p_vals)

write.table(expand.grid(pops,p_vals), 
            file = "GML_guide_pops.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)
