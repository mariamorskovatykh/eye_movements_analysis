pix2deg <- function(pix, moWidthPX,moWidthCM,distanceCM){
  # Converts a given number of pixels to degrees of visual angle
  # moWidthPX is the number of pixels that the monitor has in the horizontal
  # axis
  # moWidthCM is the width of the monitor in centimeters
  # distanceCM is the distance of the monitor to the retina
  # pix is the number of pixels you want to convert to degrees
  
  # degrees for 1 cm on the screen
  degPerCm = atan2(1,distanceCM) * 180/pi; 
  #vpixels in 1 cm on the screen
  pxPerCm = moWidthPX/moWidthCM;
  # degrees per pixel
  degPerPx = degPerCm / pxPerCm;
  # deg for given amount of pixels
  deg=degPerPx*pix;
  return(deg)
}

