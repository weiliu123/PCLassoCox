#' Plot risk scores of patients
#'
#' @param rs Risk score of patients
#' @param y.status A binary variable with 1 indicating that the event has
#'  occurred and 0 indicating (right) censoring.
#' @param cutoff The cut-off that is used to divided the patients into two
#' groups according to \code{rs}.
#' @param col.dead The color that corresponds to 1 in \code{y.status}.
#' @param col.alive The color that corresponds to 0 in \code{y.status}.
#' @param legend.position Legend position.
#' @importFrom stats median
#' @export
#'
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#'
#' cv.fit1 <- cv.PCLasso(x = GBM$GBM.train$Exp,
#'                       y = GBM$GBM.train$survData,
#'                       group = PCGroup,
#'                       nfolds = 5)
#'
#' # predict risk scores of samples in GBM.test
#' s <- predict(object = cv.fit1, x = GBM$GBM.test$Exp, type="link",
#'              lambda=cv.fit1$cv.fit$lambda.min)
#'
#' plotRS(rs = s, y.status = GBM$GBM.test$survData$status)
plotRS <- function(rs, y.status, cutoff = median(rs),
                   col.dead = "#E64B35",
                   col.alive = "#00A087",
                   legend.position = c(0.2,0.9)){
    rs.rank <- rank(rs, ties.method = "random")
    df.rs <- data.frame(x=rs.rank, y=rs,
                        status=factor(y.status,
                                      levels = c(1, 0),
                                      labels = c("Dead", "Alive")))

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    p = ggplot2::ggplot(df.rs, ggplot2::aes_string(x = "x", y = "y"))
    p + ggplot2::geom_col(ggplot2::aes_string(fill = "status"),
                 position = "dodge")+   # 直接使用数值画barplot时，用geom_col(),而不是geom_bar()
        ggplot2::geom_vline(xintercept = length(which(rs <= cutoff)), linetype = "dashed",
                   colour = "grey", size = 0.5)+
        ggplot2::labs(x="Patients (Increasing risk score)", y="Risk score")+
        ggplot2::scale_fill_manual(values = c("Alive" = col.alive, "Dead" = col.dead))+
        ggplot2::theme_classic()+
        ggplot2::theme(
            legend.position = legend.position,
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank(),
            legend.key =  ggplot2::element_blank())
}


#' Plot Kaplan-Meier curve
#'
#' @param rs Risk score of patients
#' @param y The time-to-event outcome, as a two-column matrix or \code{Surv}
#' object. The first column should be time on study (follow up time); the second
#' column should be a binary variable with 1 indicating that the event has
#' occurred and 0 indicating (right) censoring.
#' @param cutoff The cut-off that is used to divided the patients into two
#' groups according to \code{rs}.
#' @param col.low Color used to indicate low-risk patients.
#' @param col.high Color used to indicate high-risk patients.
#'
#' @importFrom stats median
#' @export
#'
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#'
#' cv.fit1 <- cv.PCLasso(x = GBM$GBM.train$Exp,
#'                       y = GBM$GBM.train$survData,
#'                       group = PCGroup,
#'                       nfolds = 5)
#'
#' # predict risk scores of samples in GBM.test
#' s <- predict(object = cv.fit1, x = GBM$GBM.test$Exp, type="link",
#'              lambda=cv.fit1$cv.fit$lambda.min)
#'
#' plotKMCurve(rs = s, y = GBM$GBM.test$survData)
#'
# survminer (>= 0.4.9)
plotKMCurve <- function(rs, y, cutoff = median(rs),
                        col.low = "#1F77B4",
                        col.high = "#FF7F0E"){
    risk <- rep(0, length = length(rs))
    risk[which(rs > cutoff)] <- 1

    subTypeData <- data.frame(risk = risk,
                              status = y$status,
                              time = y$time)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package \"survival\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("survminer", quietly = TRUE)) {
        stop("Package \"survminer\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    survminer::ggsurvplot(
        fit = survival::survfit(survival::Surv(time, status)~risk,
                                data = subTypeData),
        data = subTypeData,
        size = 1,                 # change line size
        palette =
            c(col.low, col.high),# custom color palettes
        conf.int = FALSE,          # Add confidence interval
        pval = TRUE,              # Add p-value
        risk.table = TRUE,        # Add risk table
        risk.table.col = "strata",# Risk table color by groups
        title = NULL,
        xlab = "Survival time",   # customize X axis label.
        legend.labs =
            c("Low risk", "High risk"),    # Change legend labels
        risk.table.height = 0.25, # Useful to change when you have multiple groups
        ggtheme = ggplot2::theme_classic(),      # Change ggplot2 theme
        tables.theme = ggplot2::theme_classic()
    )
}

# number: increment of the sequence.
#' Title
#'
#' @param rs Risk score of patients
#' @param y The time-to-event outcome, as a two-column matrix or \code{Surv}
#' object. The first column should be time on study (follow up time); the second
#' column should be a binary variable with 1 indicating that the event has
#' occurred and 0 indicating (right) censoring.
#' @param npoint Integer: number of time points for calculating time-dependent
#' AUC.
#' @param col Color of the time-dependent AUC curve.
#'
#' @export
#'
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#'
#' cv.fit1 <- cv.PCLasso(x = GBM$GBM.train$Exp,
#'                       y = GBM$GBM.train$survData,
#'                       group = PCGroup,
#'                       nfolds = 5)
#'
#' # predict risk scores of samples in GBM.test
#' s <- predict(object = cv.fit1, x = GBM$GBM.test$Exp, type="link",
#'              lambda=cv.fit1$cv.fit$lambda.min)
#'
#' plotTDAUC(rs = s, y = GBM$GBM.test$survData)
plotTDAUC <- function(rs, y, npoint = 50, col = "red"){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("survminer", quietly = TRUE)) {
        stop("Package \"survminer\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    timeSeq <- seq(from = 1, to=(max(y$time) + 1), by= max(y$time)/npoint)

    aucSeq <- rep(0, length=length(timeSeq))
    for(i in 1:length(timeSeq)){
        perf <- survcomp::tdrocc(x = rs,
                                surv.time = y$time,
                                surv.event = y$status,
                                time = timeSeq[i], cutpts = NA, na.rm = TRUE)
        aucSeq[i] <- perf$AUC
    }

    df.auc <- data.frame(
        time = timeSeq,
        auc = aucSeq)

    p = ggplot2::ggplot(df.auc,
                        ggplot2::aes_string(x = "time", y = "auc"))
    p + ggplot2::geom_line(color = col)+
        ggplot2::labs(x="Survival time",
                      y="Time-dependent AUC")+
        ggplot2::ylim(max(min(aucSeq)-0.3, 0),
                      min(max(aucSeq)+0.1, 1)) +
        ggplot2::theme_classic()+
        ggplot2::theme(
            legend.position = c(0.8,0.2),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank(),
            legend.key =  ggplot2::element_blank())
}
