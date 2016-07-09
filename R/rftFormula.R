#'
#'
#'
#' @export rftFormula
rftFormula <- (formula, imgData) {
  # get term specifying response/image group
  y <- all.vars(formula)[1]

  # isolate right hand side
  RHS <- delete.response(formula)

  X <- model.matrix(RHS, data = imgData$demog[imgData[y]$bool, ])
}
