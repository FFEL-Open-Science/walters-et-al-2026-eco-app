# This code takes the raw count data for each size class and formats
# them into a matrix, where the rows are survey days and the columns
# are years. All matrices are combined in an array where each slice is 
# a size class.

library(tidyverse)

count_df <- read_csv("data/raw/rb_counts.csv")

size_classes <- unique(count_df$size_class)

count_array <- array(NA, dim = c(length(unique(count_df$survey_day)),
                                 length(unique(count_df$year)),
                                 length(unique(count_df$size_class))))
                     

for (i in 1:length(size_classes)) {

  count_tmp <- filter(count_df, size_class == size_classes[i]) |>
    ungroup() |>
    mutate(year = factor(year),
           survey_day = factor(survey_day)) |>
    dplyr::select(year, survey_day, count) |>
    pivot_wider(names_from = year, values_from = count) |>
    dplyr::select(-survey_day) |>
    as.matrix()
  
  count_array[ , , i] <- count_tmp
  
}

saveRDS(count_array, "data/processed/count_data.rds")
