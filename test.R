{
   devtools::load_all()
   library(tidyverse)
   dat <- read_csv("~/Dropbox/Personal/Dissertation/data/data_rmjlca.csv")
   dat1 <- dat %>% mutate(
      ESMK_98 = ysmoke1998,
      FSMK_98 = msmoke1998 > 0,
      DSMK_98 = msmoke1998 == 30,
      HSMK_98 = pdsmoke1998 >= 10,
      EDRK_98 = ydrink1998,
      CDRK_98 = mdrink1998 > 0,
      WDRK_98 = mdrink1998 >= 5,
      BDRK_98 = bdrink1998 > 0,
      EMRJ_98 = ymari1998,
      CMRJ_98 = mmari1998 > 0,
      OMRJ_98 = mmari1998 >= 10,
      SMRJ_98 = smari1998 > 0,
      ESMK_03 = ysmoke2003,
      FSMK_03 = msmoke2003 > 0,
      DSMK_03 = msmoke2003 == 30,
      HSMK_03 = pdsmoke2003 >= 10,
      EDRK_03 = ydrink2003,
      CDRK_03 = mdrink2003 > 0,
      WDRK_03 = mdrink2003 >= 5,
      BDRK_03 = bdrink2003 > 0,
      EMRJ_03 = ymari2003,
      CMRJ_03 = mmari2003 > 0,
      OMRJ_03 = mmari2003 >= 10,
      SMRJ_03 = smari2003 > 0,
      ESMK_08 = ysmoke2008,
      FSMK_08 = msmoke2008 > 0,
      DSMK_08 = msmoke2008 == 30,
      HSMK_08 = pdsmoke2008 >= 10,
      EDRK_08 = ydrink2008,
      CDRK_08 = mdrink2008 > 0,
      WDRK_08 = mdrink2008 >= 5,
      BDRK_08 = bdrink2008 > 0,
      EMRJ_08 = ymari2008,
      CMRJ_08 = mmari2008 > 0,
      OMRJ_08 = mmari2008 >= 10,
      SMRJ_08 = smari2008 > 0
   ) %>% select(SEX = sex, RACE,
                ESMK_98, FSMK_98, DSMK_98, HSMK_98,
                EDRK_98, CDRK_98, WDRK_98, BDRK_98,
                EMRJ_98, CMRJ_98, OMRJ_98, SMRJ_98,
                ESMK_03, FSMK_03, DSMK_03, HSMK_03,
                EDRK_03, CDRK_03, WDRK_03, BDRK_03,
                EMRJ_03, CMRJ_03, OMRJ_03, SMRJ_03,
                ESMK_08, FSMK_08, DSMK_08, HSMK_08,
                EDRK_08, CDRK_08, WDRK_08, BDRK_08,
                EMRJ_08, CMRJ_08, OMRJ_08, SMRJ_08)


nlsy97 <- dat1 %>% mutate_at(vars(ends_with("98"), ends_with("03"), ends_with("08")),
                   ~factor(as.logical(.), levels = c(TRUE, FALSE), labels = c("Yes", "No"))) %>%
   mutate(SEX = factor(SEX, labels = c("Male", "Female")),
          RACE = factor(RACE, labels = c("Black", "Hispanic", "Mixed Race (Non-Hispanic)",
                                "Non-Black / Non-Hispanic")))
str(nlsy97)
save(gss1620, file = "data/gss1620.rda")
save(nlsy97, file = "data/nlsy97.rda")

   lcas <- slcm(
      smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
      smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
      smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008
   ) %>% estimate(dat1)
   lcas2 <- slcm(
      smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
      smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
      smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008,
      constraints = c("smoke1", "smoke2", "smoke3")
   ) %>% estimate(dat1)
   lcpa <- slcm(
      smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
      smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
      smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008,
      pf[3] ~ smoke1 + smoke2 + smoke3,
      constraints = c("smoke1", "smoke2", "smoke3")
   ) %>% estimate(dat1)
}
sim1 <- simulate(lcpa, 500)
fitted <- lcpa %>% estimate(data = sim1$response)
param(fitted)

gof(lcpa, test = "boot")

lca2 <- slcm(
   smoke1[2] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998
) %>% estimate(dat1)
lca3 <- slcm(
   smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998
) %>% estimate(dat1)
lca4 <- slcm(
   smoke1[4] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998
) %>% estimate(dat1)

debug(gof)
gof(lca2, lca3, lca4, test = "boot")

a <- lcas %>% regress(smoke1 ~ RACE, dat1)
debug(gof)
gof(lcpa)
class(lcpa)
lcpa <- slcm(
   smoke[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998
)
gof(object, lcas, lcpa, jlcpa)
logLik(object)
a <- poLCA::poLCA(cbind(ysmoke1998, msmoke1998, dsmoke1998, bsmoke1998) + 1 ~ 1, dat1, nclass = 3)
param(object1)
object1 <- lcpa %>% estimate(dat1)
object2 <- jlcpa %>% estimate(dat1)
objects <- list(object1, object2)

gof(object2, test = "chisq")
pp;param(object2)
gof(object, object1, object2, test = "boot")

object <- lcpa
data <- dat1
lcpa$num
# %>% fit(dat1)
param(object2)
jlcpa <- slcm(
   smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
   drink1[3] ~ ydrink1998 + mdrink1998 + ddrink1998 + bdrink1998,
   drink2[3] ~ ydrink2003 + mdrink2003 + ddrink2003 + bdrink2003,
   smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
   mari1[3] ~ ymari1998 + mmari1998 + dmari1998 + smari1998,
   mari2[3] ~ ymari2003 + mmari2003 + dmari2003 + smari2003,
   drink3[3] ~ ydrink2008 + mdrink2008 + ddrink2008 + bdrink2008,
   smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008,
   mari3[3] ~ ymari2008 + mmari2008 + dmari2008 + smari2008,
   jc1[5] ~ smoke1 + drink1 + mari1,
   jc2[5] ~ smoke2 + drink2 + mari2,
   jc3[5] ~ smoke3 + drink3 + mari3,
   pf[4] ~ jc1 + jc2 + jc3,
   constraints = list(
      c("smoke1", "smoke2", "smoke3"),
      c("drink1", "drink2", "drink3"),
      c("mari1", "mari2", "mari3"),
      c("jc1 ~ smoke1", "jc2 ~ smoke2", "jc3 ~ smoke3"),
      c("jc1 ~ drink1", "jc2 ~ drink2", "jc3 ~ drink3"),
      c("jc1 ~ mari1", "jc2 ~ mari2", "jc3 ~ mari3")
   )
) %>% estimate(dat1)
summary(jlcpa)
jlcpa <- slcm(
   smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
   drink1[3] ~ ydrink1998 + mdrink1998 + ddrink1998 + bdrink1998,
   drink2[3] ~ ydrink2003 + mdrink2003 + ddrink2003 + bdrink2003,
   smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
   mari1[3] ~ ymari1998 + mmari1998 + dmari1998 + smari1998,
   mari2[3] ~ ymari2003 + mmari2003 + dmari2003 + smari2003,
   drink3[3] ~ ydrink2008 + mdrink2008 + ddrink2008 + bdrink2008,
   smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008,
   mari3[3] ~ ymari2008 + mmari2008 + dmari2008 + smari2008,
   jc1[5] ~ smoke1 + drink1 + mari1,
   jc2[5] ~ smoke2 + drink2 + mari2,
   jc3[5] ~ smoke3 + drink3 + mari3,
   pf[4] ~ jc1 + jc2 + jc3,
   constraints = list(
      c("smoke1", "smoke2", "smoke3"),
      c("drink1", "drink2", "drink3"),
      c("mari1", "mari2", "mari3"),
      c("jc1 ~ smoke1", "jc2 ~ smoke2", "jc3 ~ smoke3"),
      c("jc1 ~ drink1", "jc2 ~ drink2", "jc3 ~ drink3"),
      c("jc1 ~ mari1", "jc2 ~ mari2", "jc3 ~ mari3")
   )
) %>% estimate(dat1)
plot(lcas, abbreviation = TRUE)
object <- jlta
summary(jlta)
jlta <- slcm(
   smoke1[3] ~ ysmoke1998 + msmoke1998 + dsmoke1998 + bsmoke1998,
   jc3[5] ~ smoke3 + drink3 + mari3,
   drink1[3] ~ ydrink1998 + mdrink1998 + ddrink1998 + bdrink1998,
   smoke3[3] ~ ysmoke2008 + msmoke2008 + dsmoke2008 + bsmoke2008,
   jc1 ~ jc2,
   smoke2[3] ~ ysmoke2003 + msmoke2003 + dsmoke2003 + bsmoke2003,
   drink2[3] ~ ydrink2003 + mdrink2003 + ddrink2003 + bdrink2003,
   jc1[5] ~ smoke1 + drink1 + mari1,
   mari2[3] ~ ymari2003 + mmari2003 + dmari2003 + smari2003,
   mari1[3] ~ ymari1998 + mmari1998 + dmari1998 + smari1998,
   jc2 ~ jc3,
   drink3[3] ~ ydrink2008 + mdrink2008 + ddrink2008 + bdrink2008,
   mari3[3] ~ ymari2008 + mmari2008 + dmari2008 + smari2008,
   jc2[5] ~ smoke2 + drink2 + mari2,
   constraints = list(
      c("jc1 ~ smoke1", "jc2 ~ smoke2", "jc3 ~ smoke3"),
      c("smoke1", "smoke2", "smoke3"),
      c("drink1", "drink2", "drink3"),
      c("mari1", "mari2", "mari3"),
      c("jc1 ~ drink1", "jc2 ~ drink2", "jc3 ~ drink3"),
      c("jc1 ~ mari1", "jc2 ~ mari2", "jc3 ~ mari3")
   )
) %>% estimate(dat1)
summary(jlta)
estimate(jlcpa, data = dat1)

ll <- -Inf
for (i in 1:20) {
   ff <- jlcpa %>% fit(data = dat1, control = list(em.iterlim = 2000))
   if (logLik(ff) > ll)
      gf <- ff
}
estimates(gf, round = 3)


jlca2$args
estimates(jlca2)

debug(fit)
lcpa %>% parameter_numbers()

debug(fit)
f <- lcpa %>% fit(restriction = paste0("rho", c(1:2)))


# , control = list(em.iterlim = 5000))
lcpa %>% parameter_numbers()
debug(fit)
f <- lcpa %>% fit(restriction = paste0("tau", c(2, 3, 17, 18, 32, 33, 4, 6, 7,8,19,21,22,23,34,36,37,38)))

