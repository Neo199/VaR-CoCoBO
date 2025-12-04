# sample_models.R

#
sample_models <- function(n_models, n_vars) {
  tmp <- sample.int(2L, size = n_models * n_vars, replace = TRUE) - 1L
  dim(tmp) <- c(n_models, n_vars)
  tmp
}

# output <- sample_models(3, 3)
# print(output)