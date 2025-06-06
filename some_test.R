
T_final <- 2

overkill_time_power <- 4
# overkill_h_power <- 10
TIME_POWERS <- overkill_time_power:0
time_steps <- 0.1 * 2^-TIME_POWERS
hs <- sqrt(time_steps) # 0.1 * 2^-c(overkill_h_power:0)

overkill_time_step <- time_steps[1]
overkill_h <- hs[1]
by_vector <- 2^rev(TIME_POWERS)

overkill_time_seq <- seq(0, T_final, by = overkill_time_step)

for (j in 1:length(time_steps)) {
  time_step <- time_steps[j]
  coarse_indices <- seq(1, length(overkill_time_seq), by = by_vector[j])
  time_seq <- overkill_time_seq[coarse_indices]
}



index_vector <- integer(length(overkill_time_seq))

for (k in 1:length(time_seq)) {
  start_idx <- which(overkill_time_seq == time_seq[k])
  end_idx <- if (k < length(time_seq)) {
    which(overkill_time_seq == time_seq[k + 1]) - 1
  } else {
    length(overkill_time_seq)
  }
  index_vector[start_idx:end_idx] <- k
}
