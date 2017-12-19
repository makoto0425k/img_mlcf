

lnglat.dist <- function(x1, y1, x2, y2, a = 6378137.0, b = 6356752.314140) { 

	deg2rad <- function(x){ 
		x * pi / 180 
	} 

	dy <- deg2rad(y1 - y2) 
	dx <- deg2rad(x1 - x2) 
	my <- deg2rad((y1 + y2) / 2) 
	e2 <- (a^2 - b^2) / a^2 
	Mnum <- a * (1 - e2) 
	W <- sqrt(1 - e2 * sin(my)^2) 
	M <- Mnum / W^3 
	N <- a / W 
	d <- sqrt((dy * M)^2 + (dx * N * cos(my))^2) 
	return(d) 
}


latlng.dist.t <- function() {
	cat(latlng.dist(139.727019, 35.626258, 	139.726492, 35.625807))
}
	
