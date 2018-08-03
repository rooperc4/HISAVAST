#' haul

#' This is a data set containing RACEBASE data extract from the GOA and AI

#' @format A data frame with 13027 rows and 33 variables:
#' \describe{
#'   \item{HAULJOIN}{unique haul identifier from RACEBASE}
#'   \item{LONGITUDE}{midpoint of haul latitude}
#'   \item{LATITUDE}{midpoint of haul longitude}
#'   \item{Year}{year of survey}
#'   \item{BOTTOM_DEPTH}{bottom depth, in meters}
#'   \item{GEAR_TEMPERATURE}{bottom water temperature, in degrees C}
#'   \item{WIRE_LENGTH}{length of wire used in the tow, in fathoms}
#'   \item{NET_WIDTH}{average net width during tow, in m}
#'   \item{DISTANCE_FISHED}{distance fished during tow in km}
#'   \item{REGION}{AI or GOA}
#'   \item{CRUISE}{cruise identifier}
#'   \item{STRATUM}{stratum designation for haul}
#'   ...
#' }
#' @source {RACEBASE}
"haul"

#' catch

#' This is a data set containing RACEBASE data extract from the GOA and AI

#' @format A data frame with 315048 rows and 10 variables:
#' \describe{
#'   \item{HAULJOIN}{unique haul identifier from RACEBASE}
#'   \item{SPECIES_CODE}{numeric species code}
#'   \item{WEIGHT}{weight of catch in survey haul}
#'   \item{NUMBER_FISH}{number of fish in the catch}
#'   \item{VESSEL}{numeric survey vessel identifier}
#'   \item{HAUL}{haul number}

#'   ...
#' }
#' @source {RACEBASE}
"catch"
