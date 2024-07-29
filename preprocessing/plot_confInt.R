Kt = 20
P = 3
var_b_given_d = pUp_results$var_b_given_d
var_b = pUp_results$var_b
spline_basis = pUp_results$spline_basis

var_tricep = spline_basis %*% as.matrix(var_b[41:60, 41:60]) %*% t(spline_basis)
var_bicep = spline_basis %*% as.matrix(var_b[21:40, 21:40]) %*% t(spline_basis)
var_Bt0 = spline_basis %*% as.matrix(var_b[1:20, 1:20]) %*% t(spline_basis)


var_tricep_givenD = spline_basis %*% as.matrix(var_b_given_d[41:60, 41:60]) %*% t(spline_basis)
var_bicep_givenD = spline_basis %*% as.matrix(var_b_given_d[21:40, 21:40]) %*% t(spline_basis)
var_Bt0_givenD = spline_basis %*% as.matrix(var_b_given_d[1:20, 1:20]) %*% t(spline_basis)


bet_df = tibble(beta0_hat = pUp_results$beta$beta0,
                bicep = pUp_results$beta$beta1,
                tricep = pUp_results$beta$beta2,
                time = pUp_results$beta$time,
                sd_beta0 = sqrt(diag(var_Bt0)),
                sd_bicep = sqrt(diag(var_bicep)),
                sd_tricep = sqrt(diag(var_tricep)),
                sd_beta0_givenD = sqrt(diag(var_Bt0_givenD)),
                sd_bicep_givenD = sqrt(diag(var_bicep_givenD)),
                sd_tricep_givenD = sqrt(diag(var_tricep_givenD)),
                ) %>%
  select(time, everything())

p0 = bet_df %>%
  select(time, contains("beta0")) %>%
  pivot_longer(sd_beta0:sd_beta0_givenD, names_to = "method", values_to = "sd") %>%
  mutate(method = ifelse(method == "sd_beta0", "sd", "sd_givenD")) %>%
  ggplot(aes(time, beta0_hat)) +
  geom_line() +
  geom_line(aes(y = beta0_hat - 2 * 1.96 * sd), linetype = 2) +
  geom_line(aes(y = beta0_hat + 2 * 1.96 * sd), linetype = 2)  +
  facet_wrap(~method) +
  ggtitle("beta0")

p1 = bet_df %>%
  select(time, contains("bicep")) %>%
  pivot_longer(sd_bicep:sd_bicep_givenD, names_to = "method", values_to = "sd") %>%
  mutate(method = ifelse(method == "sd_bicep", "sd", "sd_givenD")) %>%
  ggplot(aes(time, bicep)) +
  geom_line() +
  geom_line(aes(y = bicep - 2 * 1.96 * sd), linetype = 2) +
  geom_line(aes(y = bicep + 2 * 1.96 * sd), linetype = 2)  +
  facet_wrap(~method) +
  ggtitle("bicep")

p2 = bet_df %>%
  select(time, contains("tricep")) %>%
  pivot_longer(sd_tricep:sd_tricep_givenD, names_to = "method", values_to = "sd") %>%
  mutate(method = ifelse(method == "sd_tricep", "sd", "sd_givenD")) %>%
  ggplot(aes(time, tricep)) +
  geom_line() +
  geom_line(aes(y = tricep - 2 * 1.96 * sd), linetype = 2) +
  geom_line(aes(y = tricep + 2 * 1.96 * sd), linetype = 2)  +
  facet_wrap(~method) +
  ggtitle("tricep")
