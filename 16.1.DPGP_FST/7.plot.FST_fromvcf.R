### done in local computer
### 
library(tidyverse)
library(vroom)
library(foreach)

fst.files = list.files("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.1.DPGP_FST/fst_folder")

metadata = str_split(fst.files, "\\.")

dat.fst = foreach(i=fst.files[-4], .combine = "rbind")%do%{
  k=which(fst.files == i)
  print(i)
  tmp = vroom(paste("./fst_folder", i, sep = "/")) %>%
    mutate(W=metadata[[k]][2], S=metadata[[k]][3])
}

dat.fst %>%
  ggplot(
    aes(
      x= c(BIN_START+BIN_END)/2/1e6,
      y= WEIGHTED_FST
    )
  ) + 
  geom_vline(xintercept = 5192177/1e6, color = "blue") +
  geom_vline(xintercept = 2225744/1e6, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 13154180/1e6, color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(paste(W,S)~., ncol =1)

  