for (i in 1:nrow(poolB)){
  print(poolB[i,4])
  amp = poolB[i,4]
  if(amp %in% main[,4]){
   # print("pool2")
    main[main$Amplicon_ID==amp,5] = poolB[poolB$Amplicon_ID==amp,5]
    }
}