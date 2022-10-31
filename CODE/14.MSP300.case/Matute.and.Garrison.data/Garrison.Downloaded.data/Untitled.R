fst.2l.sim <- vroom("/Users/jcbnunez/Downloads/sim.mut.haps.fst.windowed.weir.fst")

fst.2l.sim %>%
  ggplot(aes(
    x=BIN_START,
    y=MEAN_FST
  )) + 
  geom_line() +
  ggtitle("Simulans FST - T/T vs G/G homozygous")