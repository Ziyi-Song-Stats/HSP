
##
## input curvature_index is "MAD", "AD", "AUC"
##
get.rawdata = function(curvature_index, condition){
  mousetrack = read.csv("/Users/ziyisong/Desktop/Final\ Mouse\ Data\ Application/measures_per_sub_45subs.csv")
  # delete subjects ID 77 and 166
  delete.ID = c(which(mousetrack$Subject == 77), which(mousetrack$Subject == 166))
  mousetrack = mousetrack[-delete.ID, ]
  mousetrack = data.matrix(mousetrack)
  #now we have finalized 43 subjects in total
  
  rawdata = data.frame(matrix(ncol=12, nrow=43))
  colnames(rawdata) = c("CP_NP", "CP_PH", "CP_UH",
                        "NP_CP", "NP_PH", "NP_UH",
                        "PH_CP", "PH_NP", "PH_UH",
                        "UH_CP", "UH_NP", "UH_PH")
  for (i in 1:43){
    rawdata[i, ] = mousetrack[ ,curvature_index][(12*(i-1) + 1) : (12*i)]
  }
  rawdata = data.matrix(rawdata)
  
  
  if (condition == "NP"){
    rawdata = cbind(rawdata[,"CP_NP"], rawdata[,"PH_NP"], rawdata[,"UH_NP"],
                    rawdata[,"NP_CP"], rawdata[,"NP_PH"], rawdata[,"NP_UH"])
    for (i in 1:43){
      rawdata[i, ] = ( rawdata[i, ] - mean(rawdata[i,]) ) / sd(rawdata[i,])
    }
  }
  
  if (condition == "CP"){
    rawdata = cbind(rawdata[,"CP_NP"], rawdata[,"PH_CP"], rawdata[,"UH_CP"],
                    rawdata[,"NP_CP"], rawdata[,"CP_PH"], rawdata[,"CP_UH"])
    for (i in 1:43){
      rawdata[i, ] = ( rawdata[i, ] - mean(rawdata[i,]) ) / sd(rawdata[i,])
    }
  }
  
  rawdata = t(rawdata)
  
  return(rawdata)
}



