library(thunder)

years <- seq(1980,2024)
days <- seq(1:31)

# Da 1958 al 2016: Udine Rivolto, wmo_id = 14044; dal 2016 Udine Campoformido, wmo_id = 16045

cape_udine_july <- matrix(NA, nrow=length(years), ncol=length(days))

for(i in 1:length(years)){
  for(j in 1:length(days)){
    print(c(i,j))
    if(years[i]<2016) wmo_id = 16044
    if(years[i]>=2016) wmo_id = 16045
    df <- get_sounding(wmo_id, years[i], 07, days[j], 12, metadata = FALSE)
    cape_udine_july[i,j] <- sounding_compute(df$pressure, df$altitude, df$temp, df$dpt,
                                        df$wd, df$ws, accuracy = 2, interpolate_step = 5)[1]
  }
}



