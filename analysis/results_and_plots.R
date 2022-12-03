load("blr.rdata")
begin = 25001
end = 50000

par(mfrow = c(1, 2))
hist(blr$tau_p_samples[begin:end], xlab = "tau_hat", main = "")
plot(blr$tau_p_samples[begin:end], type = "l", xlab = "estimation steps", ylab = "tau_hat samples")

mean(blr$tau_p_samples[begin:end])
quantile(blr$tau_p_samples[begin:end], c(0.025, 0.975))

par(mfrow = c(2, 4))
plot(blr$sigma_samples[begin:end, 1], type = "l", xlab = "", ylab = "samples", main = "sigma^2_0")
plot(blr$beta_0_samples[begin:end, 1], type = "l", xlab = "", ylab = "", main = "beta_{0,1}")
plot(blr$beta_0_samples[begin:end, 10], type = "l", xlab = "", ylab = "", main = "beta_{0,10}")
plot(blr$beta_0_samples[begin:end, 19], type = "l", xlab = "", ylab = "", main = "beta_{0,19}")
plot(blr$sigma_samples[begin:end, 2], type = "l", xlab = "", ylab = "samples", main = "sigma^2_1")
plot(blr$beta_1_samples[begin:end, 1], type = "l", xlab = "", ylab = "", main = "beta_{1,1}")
plot(blr$beta_1_samples[begin:end, 10], type = "l", xlab = "", ylab = "", main = "beta_{1,10}")
plot(blr$beta_1_samples[begin:end, 19], type = "l", xlab = "", ylab = "", main = "beta_{1,19}")