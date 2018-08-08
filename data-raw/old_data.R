## DATA AND MAIN FUNCTION CALL ############

#Absolute Abundance
Abs.Abundance.2005 <- data.frame(Year=c(2005), N.obs=c(6404), CV.obs=c(0.12))

#Relative Abundance
Rel.Abundance <- data.frame(Index=c(1,1,1,2,2,2),
                            Year=c(1982, 1986, 1997, 2002, 2003, 2004),
                            IA.obs=c(45,259,200, 3396, 3661, 5353),
                            CV.IA.obs=c(0.91, 0.59,0.64,0.142,0.131,0.128))

#Count data, NOT USED
Count.Data <- data.frame(Index=c(1,1,1,1,1,1,1),
                         Year=seq(2000, 2006, 1),
                         IA.obs=c(0.11, 0.15, 0.12, 0.22, 0.20, 0.14, 0.19),
                         CV.IA.obs=c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2))

#Catch Series
Catch.data <- data.frame(Year=seq(1901,2010),
                         Catch=c(0, 0, 0, 180, 288, 240, 1261, 1849, 3391,
                                 6468, 5832, 2881, 999, 1155, 1697, 447, 121,
                                 129, 111, 102, 9, 364, 133, 266, 254, 7, 0,
                                 19, 51, 107, 18, 23, 132, 57, 48, 105, 242, 0,
                                 2, 36, 13, 0, 4, 60, 238, 30, 35, 48, 83, 698,
                                 45, 34, 140, 44, 96, 167, 61, 16, 15, 27, 13,
                                 24, 12, 0, 52, 0, 189, 0, 0, 0, 0, 2, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0))

devtools::use_data(Abs.Abundance.2005,
                   Rel.Abundance,
                   Count.Data,
                   Catch.data)
