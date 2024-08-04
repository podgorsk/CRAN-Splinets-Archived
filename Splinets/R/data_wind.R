#' @title Data on wind direction and speed. 
#'
#' @description NASA/POWER CERES/MERRA2 Native Resolution Hourly Data 
#' \itemize{
#' \item Dates: 01/01/2015 through 03/05/2015 
#' \item Location: Latitude  25.7926   Longitude -80.3239 
#' \item Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = 5.4 meters
#' }
#' Data frame fields: 
#' \itemize{
#' \item \code{YEAR} -- Year of a measurement
#' \item \code{MO} -- Month of a measurement
#' \item \code{DY} -- Day of a measurement
#' \item \code{HR} -- Hour of a measurement
#' \item \code{WD10M} -- MERRA-2 Wind Direction at 10 Meters (Degrees) 
#' \item \code{WS50M} -- MERRA-2 Wind Speed at 50 Meters (m/s) 
#' \item \code{WD50M} -- MERRA-2 Wind Direction at 50 Meters (Degrees) 
#' \item \code{WS10M} -- MERRA-2 Wind Speed at 10 Meters (m/s) 
#' }
#' 
#' @docType data
#'
#' @usage data(wind)
#'
#' @format numerical \code{1536 x 8} dataframe: \code{wind} 
#'
#' @keywords datasets
#'
#' @section References:
#' The data was obtained from the National Aeronautics and Space 
#' Administration (NASA) Langley Research Center (LaRC) Prediction of 
#' Worldwide Energy Resource (POWER) Project funded through the NASA
#' Earth Science/Applied Science Program.
#' \url{https://power.larc.nasa.gov/data-access-viewer/}

#' @example R/Examples/ExWind.R
#' 
"wind"



