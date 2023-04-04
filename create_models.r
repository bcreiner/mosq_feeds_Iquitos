mod_base <- glm(cbind(num_pos, num_neg) ~ viremia + DoI + EIP, data=df_D2, family = binomial())
mod_base_re <- glmer(cbind(num_pos, num_neg) ~ viremia + DoI + EIP + (1|part_id), data=df_D2, family = binomial())
mod_factor_bin <- glm(cbind(num_pos, num_neg) ~ viremia + as.factor(DoI) + bin_EIP, data=df_D2, family = binomial())
mod_factor_bin_int <- glm(cbind(num_pos, num_neg) ~ viremia * as.factor(DoI) + bin_EIP, data=df_D2, family = binomial())

mod_gam <- gam(cbind(num_pos, num_neg) ~ s(viremia, k=3) + s(DoI, k=3) + s(EIP, k=3), data=df_D2, family = binomial())
mod_gam_bin <- gam(cbind(num_pos, num_neg) ~ s(viremia, k=3) + s(DoI, k=3) + bin_EIP, data=df_D2, family = binomial())
mod_gam_factor_bin <- gam(cbind(num_pos, num_neg) ~ s(viremia, k=3) + as.factor(DoI) + bin_EIP, data=df_D2, family = binomial())
mod_gam_linear_bin <- gam(cbind(num_pos, num_neg) ~ viremia + s(DoI, k=3) + bin_EIP, data=df_D2, family = binomial())

mod_gam_linear_factor_bin <- gam(cbind(num_pos, num_neg) ~ viremia + as.factor(DoI) + bin_EIP, data=df_D2, family = binomial())


mod_gam_eip_10p <- gam(cbind(num_pos, num_neg) ~ viremia + s(DoI, k=3) + EIP, data=df_D2[df_D2$EIP > 10,], family = binomial())
summary(mod_gam_eip_10p)
mod_gam_eip_10p <- gam(cbind(num_pos, num_neg) ~ viremia + s(DoI, k=3) + s(EIP,k=3), data=df_D2[df_D2$EIP > 10,], family = binomial())
summary(mod_gam_eip_10p)




pos <- rep(1 : length(df_D2$num_pos), df_D2$num_pos)
neg <- rep(1 : length(df_D2$num_pos), df_D2$num_neg)

df_outcome <- data.frame(rbind(df_D2[pos,],df_D2[neg,]), outcome = c(rep(1,length(pos)), rep(0, length(neg))))

mod_scam <- scam(outcome ~ s(viremia, bs="mpi",k=4) + s(DoI, k=3) + s(EIP, bs="mpi", k=4), data=df_outcome, family = binomial())

mod_scam_bin <- scam(outcome ~ s(viremia, bs="mpi",k=4) + s(DoI, k=3) + bin_EIP, data=df_outcome, family = binomial())

mod_scam_bin <- scam(outcome ~ viremia + s(DoI, k=3) + bin_EIP, data=df_outcome, family = binomial())


mod_2_lm <- lm(viremia ~ DoI, data = df_D2)
mod_2_gam <- gam(viremia ~ s(DoI, k = 3), data = df_D2)


scam(outcome ~ viremia + s(DoI, k=3) + s(EIP, bs="mpi", k=4), data=df_outcome, family = binomial())$aic

mod_A <- glm(cbind(num_pos, num_neg) ~ viremia + DoI + EIP, data=df_D2, family = binomial())
mod_B <- glm(cbind(num_pos, num_neg) ~ viremia + DoI + bin_EIP, data=df_D2, family = binomial())
mod_C <- gam(cbind(num_pos, num_neg) ~ viremia + DoI + s(EIP, k = 3), data=df_D2, family = binomial())
mod_D <- scam(outcome ~ viremia + DoI + s(EIP, bs="mpi", k=4), data=df_outcome, family = binomial())


mod_scam <- scam(outcome ~ s(viremia, bs="mpi",k=4) + s(DoI, k=3) + bin_EIP, data=df_outcome, family = binomial())
mod_scam$aic

anova(mod_A, mod_B)

