library(devtools)
load_all('../brms')
fit1 <- brm(bf(time ~ age + (age+1|disease*sex) + (1|patient), sigma ~ 1|disease*sex),
            iter = 20000, data = kidney, family = "gaussian",
            prior = set_prior("cauchy(0,2)", class = "sd"),
            marginalize = 'sex')
summary(fit1)
sink('marginalize2.stan')
stancode(fit1)
sink()