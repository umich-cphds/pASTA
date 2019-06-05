#' pASTA for multi-phenotype analysis
#'
#' description. TODO(youfei)
#'
#' @param p.values the p.value for each study.
#' @param study.sizes the sample size of each study.
#' @param cor the correlation matrix of the studies. For example, if each study
#'   is independent, cor should be the idenity matrix
#' @return a list containing the joint p value and the test statistic, which
#'   contains the optimal subset
#' @examples
#' # grab synthetic study for example
#' data("studies")
#' n.studies <- 5
#' study.sizes <- c(nrow(studies[[1]]), nrow(studies[[2]]), nrow(studies[[3]]),
#'                    nrow(studies[[4]]), nrow(studies[[5]]))
#' study.pvals <- rep(0, n.studies)
#' # Correlations of p-values among the studies.
#' # In this case the studies were generated independently so its just I
#' cor.matrix <- diag(1, n.studies)
#'# load the lrtest() function to conduct the likelihood ratio test
#'# Used just to generate the input p-values, not required in pASTA itself.
#'
#'library(lmtest)
#'
#'for(i in 1:n.studies) {
#'  # model with gene(G) by environment(E) interaction
#'  model <- glm(D ~ G + E + GbyE, data = studies[[i]], family = binomial)
#'  # model without G and GE interaction
#'  null.model <- glm(D ~ E, data = studies[[i]], family = binomial)
#'  # likelihood ratio test from the package lmtest
#'  study.pvals[i] = lmtest::lrtest(null.model, model)[2, 5]
#'}
#'
#'pasta <- pASTA(study.pvals, study.sizes, cor.matrix)
#'
#'pasta$p.pasta.joint
#'pasta$test.statistic$selected.subset
#' @references TODO(youfei)
#' @export
pASTA <- function(p.values, study.sizes, cor)
{
  statistic <- test.statistic(p.values, study.sizes)
  p.pasta.joint <- p.dlm(statistic$test.stat, study.sizes, cor)
  list(p.pasta.joint = p.pasta.joint, test.statistic = statistic)
}


p.dlm <- function(test.stat, study.sizes, cor)
{
  all.combn <- expand.grid(rep(list(0:1), length(study.sizes)))
  all.combn <- all.combn[-1,]
  sum(apply(all.combn, 1, function(v)
    integrate(cond.prob.z, test.stat, Inf, v, study.sizes, cor)$value)
  )
}

# p.values is a vector of study/trait-specific p-values
# study.size is a vector of sample sizes of the traits/studies
test.statistic <- function(p.values, study.size)
{
  # In multiple-phenotype analysis, sample sizes are usually the same for all traits
  # total number of studies
  n <- length(p.values)
  Z.stats <- -qnorm(p.values)
  all.combn <- expand.grid(rep(list(0:1), n))
  all.combn <- all.combn[-1,]
  # data frame with the last column being the Z.meta of subset S
  Z.df <- cbind(all.combn, NA)
  colnames(Z.df) <- c(paste0("Study", 1:n), "Z.S")
  rownames(Z.df) <- 1:nrow(all.combn)

  # calculate Z(S) over all possible non-empty subsets
  for (r in 1:nrow(all.combn)) {
    current.size <- study.size * all.combn[r,] # sample size of studies in a subset
    current.wt   <- sqrt(current.size / sum(current.size))
    Z.meta <- sum(current.wt * Z.stats)
    Z.df[r, n + 1] <- Z.meta
  }
  list(test.stat = max(Z.df$Z.S), Z.df = Z.df,
       selected.subset = all.combn[which.max(Z.df$Z.S),])
}

cond.prob.z <- function(z, subset, study.sizes, cor)
{
  current.size <- subset * study.sizes
  subset.size <- sum(current.size)
  # weights of studies in current subset
  current.wt <- sqrt(current.size / sum(current.size))
  # weights of studies in current subset
  n <- length(study.sizes)
  p <- 0
  # integrate over all studies not included in subset
  for (k in 1:n) {
    e <- rep(0, n)
    e[k] <- 1
    A <- rbind(e, current.wt)
    sigma <- A %*% cor %*% t(A)
    cond.mean <- sigma[1, 2] / sigma[2, 2] * z
    cond.sigma  <- sqrt(sigma[1, 1] - sigma[1, 2] / sigma[2, 2] * sigma[2, 1])
    a <- subset.size / (subset.size + study.sizes[k])
    if (k %in% which(subset == 0)) {
      # conditional probability given Z_gamma = z
      a <- subset.size / (subset.size + study.sizes[k])
      p <- p + log(pnorm(z * (1 - sqrt(a)) / sqrt(1 - a), cond.mean, cond.sigma))
    }
    # k-th study included in subset
    else {
      b <- subset.size / (subset.size - study.sizes[k])
      if (sum(subset) > 1)
        p <- p + log(pnorm((z * (sqrt(b) - 1)) / sqrt(b - 1), cond.mean,
                       cond.sigma, lower.tail = F))
    }
  }

  exp(p) * dnorm(z, sd = sqrt(t(current.wt) %*% cor %*% current.wt))
}
