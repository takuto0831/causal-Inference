library(CausalImpact)

set.seed(1)
x1 <- 100 + arima.sim(model = list(ar = 0.999), n = 100)
x1_fix <- scale(x1)

y <- 1.7 * x1 + rnorm(100)

y[71:100] <- y[71:100] + 10
data <- cbind(y, x1)

pre.period <- c(1, 70)
post.period <- c(71, 100)

impact <- CausalImpact(data, 
                       pre.period, 
                       post.period, 
                       model.args = list(niter=63))
plot(impact)

## get result df 
df <- impact[[1]]
df |> data.frame() |> View()

## how to get model parameter
list_ <- impact[[4]]
summary(impact)

multiple.func <- function(x) {
  c(mean = mean(x), median= median(x))
}

### result dfの point.predと同等? 
list_$posterior.samples |> 
  data.frame() |>
  sapply(FUN=multiple.func) 

## bsts model
bsts_model <- list_$bsts.model

# burn in の期間, その部分をカットして trend, regressionの係数 (Z_t)
burn <- SuggestBurn(0.1, bsts_model)
bsts_model$state.contributions[-(1:burn), ,]
 
# trend (u_t)は良い感じに求められている, regressionはbeta*x (標準化済み)を表しているはず. 
tmp2 <- data_frame(
  x = 1:100,
  trend = bsts_model$state.contributions[-(1:burn),"trend" ,] |> colMeans(),
  regression = bsts_model$state.contributions[-(1:burn),"regression" ,] |> colMeans(),
  x1 = x1, 
  x1_fix = x1_fix) |> 
  mutate(beta = regression/ x1_fix,
         maybe_y_samples = trend + regression)

tmp2 |> 
  gather("component", "value", trend, regression, x1_fix, x1) %>%
  ggplot(aes(x = x, y= value)) + 
  geom_line() + theme_bw() + 
  theme(legend.title = element_blank()) + ylab("") + xlab("") +
  facet_grid(component ~ ., scales="free") + guides(colour=FALSE) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

plot(bsts_model, "components")
bsts_model$state.specification

### bstsと同じ定義なら, ここが係数 -> 不要なら0になるはず.
# y,x1を標準化した上での観測方程式の係数行列となる.
bsts_model$coefficients
bsts_model$sigma.obs # 観測誤差

plot(bsts_model, "coefficients")
plot(bsts_model, "predictors")

# カルマンフィルターによる1期間先の予測の誤差 -> これ使えそう.
tmp <- bsts_model$one.step.prediction.errors |>
  data.frame() |> 
  sapply(FUN=mean) |> 
  as.vector()
plot(tmp)
#df |> 
#  data.frame() |> 
#  mutate(diff_hand = response-point.pred)

# こういうのも動く
summary.bsts(bsts_model)

# 関数見ながら試す.
library(assertthat)
model.args = NULL
bsts.model = NULL
post.period.response = NULL
alpha=0.05

### analysis: CausalImpact function
# https://github.com/google/CausalImpact/blob/master/R/impact_analysis.R
checked <- FormatInputForCausalImpact(data, pre.period, post.period,
                                      model.args, bsts.model,
                                      post.period.response, alpha)
data <- checked$data
pre.period <- checked$pre.period
post.period <- checked$post.period
model.args <- checked$model.args
bsts.model <- checked$bsts.model
post.period.response <- checked$post.period.response
alpha <- checked$alpha

### analysis: RunWithData function
# time points.
times <- time(data)
time(data) <- seq_len(nrow(data))

# Zoom in on data in modeling range, remove original time indices.
pre.period[1] <- max(pre.period[1], which.max(!is.na(data[, 1])))
data.modeling <- window(data, start = pre.period[1])
times.modeling <- window(times, start = pre.period[1])
if (is.null(ncol(data.modeling))) {
  dim(data.modeling) <- c(length(data.modeling), 1)
}

# Standardize all variables?
UnStandardize <- identity
if (model.args$standardize.data) {
  fit.range <- c(1, diff(pre.period) + 1)
  sd.results <- StandardizeAllVariables(data.modeling, fit.range)
  data.modeling <- sd.results$data
  UnStandardize <- sd.results$UnStandardize
}

### impact_inference: CompilePosteriorInferences
# https://github.com/google/CausalImpact/blob/master/R/impact_inference.R

y.samples <- ComputeResponseTrajectories(bsts_model)
state.samples <- GetPosteriorStateSamples(bsts_model)
point.pred <- ComputePointPredictions(y.samples, state.samples, alpha=0.05)

y.samples |> colMeans()
# trend + regressionの数値とほぼ同じになるはず.
# -> 実際のy samplesは平均じゃないから一致しない? code読めばわかる.


# Undo standardization (if any)
y.samples_ <- UnStandardize(y.samples)
point.pred_ <- UnStandardize(point.pred)
y.model_ <- UnStandardize(bsts.model$original.series)

