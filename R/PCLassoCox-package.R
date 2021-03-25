#' @title \packageTitle{PCLassoCox}
#' @description \packageDescription{PCLassoCox}
#' @details
#' \packageIndices{PCLasso} The PCLasso model accepts a gene expression matrix,
#' survival data, and protein complexes for the PCLasso model, and makes
#' predictions for new samples and identifies risk protein complexes.
#'
#' \code{PCLasso} constructs a \code{PCLasso} model based on a gene expression
#' matrix, survival data, and protein complexes.
#'
#' \code{predict.PCLasso} makes predictions from a \code{PCLasso} model.
#'
#' \code{cv.PCLasso} performs k-fold cross validations for the \code{PCLasso}
#' model with grouped covariates over a grid of values for the regularization
#' parameter \code{lambda}, and returns an optimal value for \code{lambda}.
#'
#' \code{predict.cv.PCLasso} returns predictions from a fitted \code{cv.PCLasso}
#' object, using the optimal value chosen for \code{lambda}.
#'
#' \code{plot.PCLasso} produces a plot of the coefficient paths for a fitted
#' \code{PCLasso} object.
#'
#' \code{plot.cv.PCLasso} plots the cross-validation curve from a
#' \code{cv.PCLasso} object, along with standard error bars.
#'
#'
#' @references
#' PCLasso: a protein complex-based group lasso-Cox model for accurate prognosis
#' and risk protein complex discovery. To be published.
#'
#' Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
#' lasso for integrative multi-omics analysis. Journal of computational biology:
#' a journal of computational molecular cell biology, 22, 73-84.
#'
#' @docType package
#' @name PCLassoCox
NULL
#> NULL

