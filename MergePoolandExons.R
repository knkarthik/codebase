#main is the dataframe derived from the 384well design plate from Life Tech and has the information for all the amplicons in all the pools
#In addition, it has a column with the exon numbers covered by each amplicon
#poolB is a subset of main that has only the information for well, poolB

for (i in 1:nrow(poolB)){
  #get the ampliconID
  amp = poolB[i,4]
  if(amp %in% main[,4]){
   # get the exon number corresponding to the ampliconID from poolB and add it to the corresponding exon column of main
    main[main$Amplicon_ID==amp,5] = poolB[poolB$Amplicon_ID==amp,5]
    }
}
