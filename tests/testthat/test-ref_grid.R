context("Reference grids")

pigs.lm = lm(log(conc) ~ source + factor(percent), data = pigs)
rg = ref_grid(pigs.lm)
rg1 = ref_grid(pigs.lm, at = list(source = "soy", percent = 12))

pigs.lm2 = update(pigs.lm, conc ~ source + percent, data = pigs)
rg2 = ref_grid(pigs.lm2)
rg2a = ref_grid(pigs.lm2, at = list(source = c("fish", "soy"), percent = 10))
rg2c = ref_grid(pigs.lm2, cov.reduce = FALSE)
rg2m = ref_grid(pigs.lm2, cov.reduce = min)

pigs.lm3 = update(pigs.lm2, . ~ source + source:factor(percent))
pigs = transform(pigs, sp = interaction(source, percent))
pigs.lm4 = update(pigs.lm2, . ~ source + sp)

test_that("Reference grid is constructed correctly", {
    expect_equal(nrow(rg@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(nrow(rg2@grid), 3)
    expect_equal(nrow(rg2a@grid), 2)
    expect_equal(nrow(rg2c@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(length(rg@levels), 2)
    expect_equal(rg2@levels$percent, mean(pigs$percent))
    expect_equal(rg2m@levels$percent, min(pigs$percent))
})

test_that("Bootstrap covariance is available for model-based reference grids", {
    boot.emm = emmeans(pigs.lm2, "source", bootstrap = TRUE,
                       bootstrap.n = 40, bootstrap.seed = 27)
    expect_equal(dim(boot.emm@V), c(length(boot.emm@bhat), length(boot.emm@bhat)))
    expect_equal(boot.emm@misc$bootstrap$method, "simulate/refit")
    expect_equal(boot.emm@misc$bootstrap$n, 40)
    boot.default = summary(boot.emm)
    expect_true("SE" %in% names(boot.default))
    expect_false(any(c("df", "t.ratio", "z.ratio") %in% names(boot.default)))
    boot.sum = summary(boot.emm, infer = c(FALSE, TRUE), adjust = "none")
    boot.est = boot.emm@linfct %*% t(boot.emm@misc$bootstrap$draws)
    expect_equal(boot.default$emmean, rowMeans(boot.est))
    expect_equal(boot.default$SE, apply(boot.est, 1, sd))
    expected.p = pmin(1, 2 * pmin(rowMeans(boot.est <= 0),
                                  rowMeans(boot.est >= 0)))
    expect_equal(boot.sum$p.value, expected.p)
    boot.ci = summary(boot.emm, infer = c(TRUE, FALSE), adjust = "none")
    expected.ci = t(apply(boot.est, 1, quantile,
                          probs = c(.025, .975), names = FALSE))
    expect_equal(boot.ci$lower.CL, expected.ci[, 1])
    expect_equal(boot.ci$upper.CL, expected.ci[, 2])
    boot.left = summary(boot.emm, infer = c(TRUE, FALSE), side = "<")
    expect_true(all(is.infinite(boot.left$lower.CL) & boot.left$lower.CL < 0))
    expected.left = t(apply(boot.est, 1, quantile,
                            probs = c(0, .95), names = FALSE))
    expect_equal(boot.left$upper.CL, expected.left[, 2])
    boot.right = summary(boot.emm, infer = c(TRUE, FALSE), side = ">")
    expected.right = t(apply(boot.est, 1, quantile,
                             probs = c(.05, 1), names = FALSE))
    expect_equal(boot.right$lower.CL, expected.right[, 1])
    expect_true(all(is.infinite(boot.right$upper.CL) & boot.right$upper.CL > 0))
    expect_equal(unname(bootstrap_samples(boot.emm)), unname(t(boot.est)))
    expect_equal(bootstrap_samples(boot.emm, "coefficients"),
                 boot.emm@misc$bootstrap$draws)
})

test_that("Mixed-model bootstrap methods are selectable", {
    skip_if_not_installed("lme4")
    fit = lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
    boot.emm = emmeans(fit, "Days", bootstrap = TRUE,
                       bootstrap.type = "parametric", bootstrap.n = 10,
                       bootstrap.seed = 27)
    expect_equal(boot.emm@misc$bootstrap$method, "simulate/refit")
    expect_equal(boot.emm@misc$bootstrap$n, 10)
    expect_equal(nrow(bootstrap_samples(boot.emm, "coefficients")), 10)
    boot.summary = summary(boot.emm)
    expect_true("SE" %in% names(boot.summary))
    expect_false(any(c("df", "t.ratio", "z.ratio") %in% names(boot.summary)))
    expect_error(emmeans(fit, "Days", bootstrap = TRUE,
                          bootstrap.type = "unknown"), "should be one of")
})

test_that("GEE bootstrap uses coefficient draws", {
    skip_if_not_installed("geepack")
    gee.data = data.frame(y = c(1, 0, 2, 1, 3, 1, 2, 0),
                          trt = factor(rep(c("a", "b"), 4)),
                          id = rep(1:4, each = 2))
    gee.mod = geepack::geeglm(y ~ trt, id = id, data = gee.data,
                              family = poisson("log"))
    gee.rg = ref_grid(gee.mod, bootstrap = TRUE, bootstrap.n = 40,
                      bootstrap.seed = 28)
    expect_equal(gee.rg@misc$bootstrap$method, "coefficient draws")
    expect_true(all(is.finite(gee.rg@V)))
})

test_that("GLM bootstrap resampling methods are available", {
    gaussian.glm = glm(conc ~ source, data = pigs, family = gaussian())
    residual.rg = ref_grid(gaussian.glm, bootstrap = TRUE,
                           bootstrap.type = "residual", bootstrap.n = 10,
                           bootstrap.seed = 29)
    expect_equal(residual.rg@misc$bootstrap$method, "residual/refit")
    binary.glm = glm(I(conc > median(conc)) ~ source, data = pigs,
                     family = binomial())
    case.rg = ref_grid(binary.glm, bootstrap = TRUE,
                       bootstrap.type = "case", bootstrap.n = 10,
                       bootstrap.seed = 30)
    expect_equal(case.rg@misc$bootstrap$method, "case/refit")
})

test_that("GEE case bootstrap is available", {
    skip_if_not_installed("geepack")
    gee.data = data.frame(y = c(1, 0, 2, 1, 3, 1, 2, 0),
                          trt = factor(rep(c("a", "b"), 4)),
                          id = rep(1:4, each = 2))
    gee.mod = geepack::geeglm(y ~ trt, id = id, data = gee.data,
                              family = poisson("log"))
    gee.rg = ref_grid(gee.mod, bootstrap = TRUE, bootstrap.type = "case",
                      bootstrap.n = 10, bootstrap.seed = 31)
    expect_equal(gee.rg@misc$bootstrap$method, "case/refit")
})

test_that("Gaussian GEE residual bootstrap is available", {
    skip_if_not_installed("geepack")
    gee.data = data.frame(y = c(1, 0, 2, 1, 3, 1, 2, 0),
                          trt = factor(rep(c("a", "b"), 4)),
                          id = rep(1:4, each = 2))
    gee.mod = geepack::geeglm(y ~ trt, id = id, data = gee.data,
                              family = gaussian())
    gee.rg = ref_grid(gee.mod, bootstrap = TRUE,
                      bootstrap.type = "residual", bootstrap.n = 10,
                      bootstrap.seed = 32)
    expect_equal(gee.rg@misc$bootstrap$method, "residual/refit")
})

test_that("Reference grid extras are detected", {
    expect_equal(rg@misc$tran, "log")
    expect_true(is.null(rg2@misc$tran))
    expect_true(is.null(rg2@model.info$nesting))
    expect_is(ref_grid(pigs.lm3)@model.info$nesting, "list") # see note above
    expect_is(ref_grid(pigs.lm4)@model.info$nesting, "list") # see note above
})

colnames(ToothGrowth) <- c('len', 'choice of supplement', 'dose')
model <- stats::aov(`len` ~ `choice of supplement`, ToothGrowth)

test_that("Reference grid handles variables with spaces", {
    expect_output(str(ref_grid(model, ~`choice of supplement`)), "choice of supplement")
})

# models outside of data.frames
x = 1:10
y = rnorm(10)
mod1 = with(pigs, lm(log(conc) ~ source + factor(percent)))
test_that("ref_grid works with no data or subset", {
    expect_silent(ref_grid(lm(y ~ x)))
    expect_silent(ref_grid(mod1))
})

# Multivariate models
MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
MOats.rg <- ref_grid (MOats.lm, mult.levs = list(
    trt = LETTERS[1:2], dose = as.character(1:2))
)
test_that("We can construct multivariate reference grid", {
    expect_equal(nrow(MOats.rg@grid), 72)
    expect_equal(length(MOats.rg@levels), 4)
})
MOats.nrg <- ref_grid (MOats.lm, mult.levs = list(nitro = c(0,.2,.4,.6)),
                       at = list(nitro = c(0.0001, .39998, .5)))
test_that("Fuzzy matching of numerical mult.levels works", {
    expect_equal(length(MOats.nrg@levels$nitro), 2)
    expect_equal(MOats.nrg@levels$nitro, c(0, .4), 0.001)
})

### Nuisance factors
MOats.rgn = ref_grid(MOats.lm, nuisance = "Block")
MOats.emm = emmeans(MOats.lm, ~ Variety * rep.meas)
test_that("We get same predictions with and without nuisance specs", {
    expect_equal(predict(MOats.rgn), predict(MOats.emm), 0.001)
})

# Missing and NA levels
miss.df = data.frame(x = factor(c("a", "a", "b", NA), levels = c("a", "b", "c", NA)), y = 1:4)
miss.lm = lm(y ~ x, data = miss.df)
miss.rg1 = ref_grid(miss.lm)
miss.rg2 = ref_grid(miss.lm, data = miss.df)
# Now try allowing NA levels
miss.dfa = transform(miss.df, x = factor(x, exclude = NULL))
miss.lma = lm(y ~ x, data = miss.dfa)
miss.rg1a = ref_grid(miss.lma)
miss.rg2a = ref_grid(miss.lma, data = miss.dfa)
test_that("Reference grid handles missing values", {
    expect_equal(length(miss.rg1@levels$x), 2)
    expect_equal(length(miss.rg2@levels$x), 2)
    expect_equal(length(miss.rg1a@levels$x), 3)
    expect_equal(length(miss.rg2a@levels$x), 3)
    expect_equal(inherits(with_emm_options(allow.na.levs = FALSE, ref_grid(miss.lma)), 
                 "try-error"), TRUE)
})

