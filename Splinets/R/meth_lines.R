#' @title Adding graphs of splines to a plot 
#'
#' @description  A standard method of adding splines to an existing plot. 
#' @param x \code{Splinets} object;
#' @param sID vector, specifying indices of splines in the splinet object to be plotted;
#' @param ... other standard graphical parameters;
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{plot,Splinets-method}} for graphical visualization of splines;
#' \code{\link{evspline}} for evaluation of a \code{Splinet}-object;
#'
#' @example R/Examples/ExLines.R
#' 
#' 

# setGeneric("lines.Splinets",
#            function(object, sID = NULL, ...){
#              standardGeneric("lines.Splinets")
#            }
# )


#In the most recent version the programs did not run if one used 'lines' instead of 'lines.Splinets'.
#Thus the following is used (similarily as in plot())
setMethod(
  "lines","Splinets",
  function(x, sID = NULL, ...){
    lines.spline(x, sID = sID, ...)
}
)
