#' The GBM datasets
#'
#' A dataset containing the GBM expression data and survival data.
#'
#'
#' @format A list with four GBM datasets:
#' \describe{
#'   \item{GBM.train}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{GBM.train}.
#'   \code{survData}: Survival data of \code{GBM.train}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{GBM.test}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{GBM.test}.
#'   \code{survData}: Survival data of \code{GBM.test}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{GBM.valid1}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{GBM.valid1}.
#'   \code{survData}: Survival data of \code{GBM.valid1}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{GBM.valid2}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{GBM.valid2}.
#'   \code{survData}: Survival data of \code{GBM.valid2}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'   }
#'
#' @usage data(GBM)
#' @examples data(GBM)
"GBM"

#' The LUAD datasets
#'
#' A dataset containing the LUAD expression data and survival data.
#'
#'
#' @format A list with five LUAD datasets:
#' \describe{
#'   \item{LUAD.train}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{LUAD.train}.
#'   \code{survData}: Survival data of \code{LUAD.train}. The first column is
#'   the time on study (follow up time); the second column is a binary variable
#'   with 1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{LUAD.test}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{LUAD.test}.
#'   \code{survData}: Survival data of \code{LUAD.test}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{LUAD.valid1}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{LUAD.valid1}.
#'   \code{survData}: Survival data of \code{LUAD.valid1}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{LUAD.valid2}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{LUAD.valid2}.
#'   \code{survData}: Survival data of \code{LUAD.valid2}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{LUAD.valid3}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{LUAD.valid3}.
#'   \code{survData}: Survival data of \code{LUAD.valid3}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'   }
#'
#' @usage data(LUAD)
#' @examples data(LUAD)
"LUAD"

#' The OV datasets
#'
#' A dataset containing the OV expression data and survival data.
#'
#'
#' @format A list with four OV datasets:
#' \describe{
#'   \item{OV.train}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{OV.train}.
#'   \code{survData}: Survival data of \code{OV.train}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{OV.test}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{OV.test}.
#'   \code{survData}: Survival data of \code{OV.test}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{OV.valid1}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{OV.valid1}.
#'   \code{survData}: Survival data of \code{OV.valid1}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'
#'   \item{OV.valid2}{
#'   A list with two items.
#'   \code{Exp}: Expression data of \code{OV.valid2}.
#'   \code{survData}: Survival data of \code{OV.valid2}. The first column is the
#'   time on study (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating(right)
#'   censoring.}
#'   }
#'
#' @usage data(OV)
#' @examples data(OV)
"OV"

#' Protein complexes for "PCLasso"/"cv.PCLasso"
#'
#' A list of protein complexes. The genes in each protein complex are
#' represented by EntrezID, which are consistent with the gene names in the
#' expression data in \code{GBM}, \code{LUAD}, and \code{OV}.
#'
#'
#' @format A list of 2417 protein complexes.
#' @usage data(PCGroup)
#' @examples data(PCGroup)
"PCGroup"
