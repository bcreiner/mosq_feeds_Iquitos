if (reprep){
  DF_names <- c("part_id", "pathogen", "sex", "age", "index_contact",
                "symptomatic", "febrile", "diagnosis", "sample_code",
                "collection_date", "DoI", "day_symptom", "result_deng",
                "viremia", "num_tested", "num_pos", "EIP")
  
  feed_raw <- read.csv("data/p01_feed_sep_2019.csv", as.is = TRUE)[,-1]
  
  infectiousness <- read.csv("data/infectiousness_cleaned.csv", as.is = TRUE)[,-1]
  UID <- unique(infectiousness$participant_id)
  infectiousness <- infectiousness[which(infectiousness$salivated == TRUE & !is.na(infectiousness$num_innoc) & infectiousness$diagnosis != "negative"),]
  infectiousness$collection_date <- as.Date(infectiousness$collection_date)
  infectiousness$date_fate <- as.Date(infectiousness$date_fate)
  
  infectiousness$test_duration <- infectiousness$date_fate - infectiousness$collection_date
  
  df <- data.frame(matrix(ncol = length(DF_names), nrow = 0, dimnames=list(NULL, DF_names)))
  
  for (i in 1:length(UID)){
    tmplocs <- which(infectiousness$participant_id == UID[i])
    tmpsamp <- unique(infectiousness$sample_code[tmplocs])
    for (j in 1:length(tmpsamp)){
      tmplocs2 <- which(infectiousness$sample_code == tmpsamp[j])
      tmpvals <- unique(infectiousness$test_duration[tmplocs2])
      tmpmat <- infectiousness[tmplocs2,]
      for (k in 1:length(tmpvals)){
        tmplocs3 <- which(tmpmat$test_duration == tmpvals[k])
        tmpmat2 <- tmpmat[tmplocs3,]
        wastest <- which(!is.na(tmpmat2$result_pool))
        if (length(wastest)){
          tmptest <- length(wastest)
          tmpmat3 <- tmpmat2[wastest,]
          tmppos <- length(which(tmpmat3$result != "NEG" ))
          ##
          tmp_df <- data.frame(matrix(ncol = length(DF_names), nrow = 1, dimnames=list(NULL, DF_names)))
          tmp_df$part_id <- i
          tmp_df$pathogen <- tmpmat2$pathogen[1]
          tmp_df$sex <- tmpmat2$sex[1]
          tmp_df$age <- tmpmat2$age[1]
          tmp_df$index_contact <- tmpmat2$index_contact[1]
          tmp_df$symptomatic <- tmpmat2$symptomatic[1]
          tmp_df$febrile <- tmpmat2$febrile[1]
          tmp_df$diagnosis <- tmpmat2$diagnosis[1]
          tmp_df$sample_code <- tmpmat2$sample_code[1]
          tmp_df$collection_date <- tmpmat2$collection_date[1]
          tmp_df$DoI <- tmpmat2$day_fever[1]
          tmp_df$day_symptom <- tmpmat2$day_symptom[1]
          tmp_df$result_deng <- tmpmat2$result_deng[1]
          tmp_df$viremia <- tmpmat2$engorged_titre_mean_all[1]
          tmp_df$num_tested <- tmptest
          tmp_df$num_pos <- tmppos
          tmp_df$EIP <- as.numeric(tmpmat2$test_duration[1])
          #
          df <- rbind(df, tmp_df)
        }
      }
    }
  }
  
  df$fracPos <- df$num_pos / df$num_tested
  df$num_neg <- df$num_tested - df$num_pos
  
  df <- df[which(df$pathogen != "zika"), ]
  df <- df[,c("part_id", "diagnosis", "DoI", "result_deng", 
              "viremia", "num_tested", "num_pos", "EIP", "fracPos", "num_neg")]
  
  write.csv(df, "input/manuscript_data.csv")
} else {
  df <- read.csv("input/manuscript_data.csv", as.is = TRUE)[,-1]
}

UID <- unique(infectiousness$participant_id)
if (Include_Negative_NA){
  df_D2 <- df[-which(df$result_deng == "DENV-3"),]
} else {
  df_D2 <- df[which(df$result_deng == "DENV-2"),]
}

df_D2$bin_EIP <- as.factor(ifelse(df_D2$EIP > 10, 2,1))

df_D2 <- df_D2[which(!is.na(df_D2$viremia)), ]
df_D2$part_id <- as.factor(df_D2$part_id)
