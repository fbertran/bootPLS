#' Graphical assessment of the stability of selected variables
#' 
#' This function is based on the \code{\link[bipartite]{visweb}} function from
#' the bipartite package.
#' 
#' @param matbin Matrix with 0 or 1 entries. Each row per predictor and a
#' column for every model. 0 means the predictor is not significant in the
#' model and 1 that, on the contrary, it is significant.
#' @param pred.lablength Maximum length of the predictors labels. Defaults to
#' full label length.
#' @param labsize Size of the predictors labels.
#' @param plotsize Global size of the graph.
#' @return A plot window.
#' @author Bernd Gruber with minor modifications from
#' Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See Also \code{\link{visweb}}
#' @references Vazquez, P.D., Chacoff, N.,P. and Cagnolo, L. (2009) Evaluating
#' multiple determinants of the structure of plant-animal mutualistic networks.
#' \emph{Ecology}, 90:2039-2046.
#' @keywords hplot
#' @export signpred2
#' @examples
#' set.seed(314)
#' simbin <- matrix(rbinom(200,3,.2),nrow=20,ncol=10)
#' signpred2(simbin)
signpred2 <- function(matbin, pred.lablength = max(sapply(rownames(matbin),
                     nchar)), labsize = 1, plotsize = 12)
{
  lll = 2
  text = "no"
  ncol <- ncol(matbin)
  nrow <- nrow(matbin)
  plotsize = plotsize/2.54
  mcol = max(matbin)
  if (ncol > nrow) {
    wx <- plotsize
    wy <- (plotsize)/ncol * nrow
  }
  else {
    wy <- plotsize
    wx <- (plotsize)/nrow * ncol
  }
  m.colsize = max(strwidth(colnames(matbin), units = "inches"))
  m.rowsize = max(strwidth(rownames(matbin), units = "inches"))
  cellsize = wx/ncol
  if (substr(text, 1, 1) == "i")
    s <- as.character(max(matbin))
  else s = "A"
  lettersize = strwidth(s, units = "inches")
  clratio = cellsize/lettersize
  mm <- max(m.colsize, m.rowsize)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(las = 3, mar = c(2, 2, 1, 1) + 0.1, mgp = c(2, 1, 0))
  bipartite::visweb(t(matbin), type = "None", labsize = labsize, 
                    square = "defined", def.col = c("white","grey25","green"), 
                    prednames = TRUE, clear = FALSE)
}
