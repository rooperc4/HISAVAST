#' Net position
#'
#' This function calculates the start and end net positions on the seafloor of a survey
#' tow using the wire out, depth and start and end latitude and longitude. Uses geometry
#' and direction of the tow to estimate the offset behind the vessel.
#' @param START_LAT start latitude of a survey tow
#' @param END_LAT end latitude of a survey tow
#' @param START_LON start longitude of a survey tow
#' @param END_LON end longitude of a survey tow
#' @param WIRE_OUT wire out during a survey tow 
#' @param BOTTOM_DEPTH mean bottom depth fo a survey tow 
#' @keywords net positions
#' @export
#' @examples
#' NetPosition()
NetPosition<-function(START_LAT,END_LAT,START_LON,END_LON,WIRE_OUT,BOTTOM_DEPTH){
  START_LAT<-START_LAT*pi/180
  END_LAT<-END_LAT*pi/180
  START_LON<-START_LON*pi/180
  END_LON<-END_LON*pi/180
  
  Bearing<-((atan2(sin(END_LON-START_LON)*cos(END_LAT),cos(START_LAT)*sin(END_LAT)-sin(START_LAT)*cos(END_LAT)*cos(END_LON-START_LON))*180/pi)+360)%%360
  Distance_behind<-sqrt(WIRE_OUT^2-BOTTOM_DEPTH^2)*-1
  
  START_LATC<-asin(sin(START_LAT)*cos(Distance_behind/6367449)+cos(START_LAT)*sin(Distance_behind/6367449)*cos(Bearing*pi/180))*180/pi
  END_LATC<-asin(sin(END_LAT)*cos(Distance_behind/6367449)+cos(END_LAT)*sin(Distance_behind/6367449)*cos(Bearing*pi/180))*180/pi
  
  START_LONC<-(START_LON+atan2(sin(Bearing*pi/180)*sin(Distance_behind/6367449)*cos(START_LAT),cos(Distance_behind/6367449)-sin(START_LAT)*sin(START_LATC*pi/180)))*180/pi
  END_LONC<-(END_LON+atan2(sin(Bearing*pi/180)*sin(Distance_behind/6367449)*cos(END_LAT),cos(Distance_behind/6367449)-sin(END_LAT)*sin(END_LATC*pi/180)))*180/pi
  
  return(data.frame(START_LATC,END_LATC,START_LONC,END_LONC))}