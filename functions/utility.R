#helper functions for data formatting
FormatDiffDF <- function(df, idvar, timevar, Yvar) {
  df <- df %>%
    rename(Y = {{Yvar}}, Patient = {{idvar}}, t = {{timevar}}) %>%
    arrange(Patient, t) %>%
    group_by(Patient) %>%
    mutate(V = Y - lag(Y),
           tau = t - lag(t),
           T_max = max(t),
           n_visits = n())
  
  df <- df %>% 
    dplyr::filter(!is.na(V) & !is.na(tau)) %>% 
    dplyr::select(Patient, V, tau, Y, t, T_max, n_visits) %>%
    dplyr::filter(n_visits > 1, tau>0)
  
  return(df)
}
PadList <- function(l, ...) {
  max_length <- max(sapply(l, length))
  padded_list <- lapply(l, function(x) c(x, rep(NA, max_length - length(x))))
  result_matrix <- do.call(cbind, padded_list)
  return(result_matrix)
}
MakeLong <- function(data, time, ...) {
  
  data_df <- as.data.frame(data)
  time_df <- as.data.frame(time)
  data_df$visit_num <- 1:nrow(data_df)
  time_df$visit_num <- 1:nrow(time_df)
  
  data_long <- data_df %>%
    pivot_longer(-visit_num, names_to = "Patient", values_to = "Y")
  time_long <- time_df %>%
    pivot_longer(-visit_num, names_to = "Patient", values_to = "t")
  long_tbl <- left_join(data_long, time_long, by = c("visit_num", "Patient"))
  long_tbl <- long_tbl %>% drop_na(Y, t) %>% arrange(Patient)
  
  return(long_tbl)
}
AddDiff <- function(df, Y, Patient, t, ...) {
  df %>%
    rename(Y = {{Y}}, Patient = {{Patient}}, t = {{t}}) %>%
    arrange(Patient, t) %>%
    group_by(Patient) %>%
    mutate(V = Y - lag(Y),
           tau = t - lag(t),
           T_max = max(t),
           n = n())
}