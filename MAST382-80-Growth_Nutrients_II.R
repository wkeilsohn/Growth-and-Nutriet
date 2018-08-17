
#####In this homework assignment, you will learn how to read in data from the web, 
#####and calculate underwater light fields using the Lambert-Beer Law
#######and calculate underwater fields of primary production.
######Then we will calculate growth from primary production

require(fields) ####loading the fields library (necessary for some of the plotting)
################Reading in Data

sun_light <- read.csv("http://modata.ceoe.udel.edu/public/five_year_light_data_627.csv", header=T) ###NOTE: This file can also be downloaded and loaded from your own computer

sun_light = sun_light[1:8760,]

sun_time <- paste(paste(sun_light$Year, sun_light$Month, sun_light$Day, sep = "-"), paste(sun_light$Time, "00:00", sep = ":"))
sun_time <- strptime(sun_time, format="%Y-%m-%d %H:%M:%S", tz = "GMT")
sun_time <- as.POSIXct(sun_time, origin = "1970-01-01", tz = "GMT")

PAR <- sun_light[,5]*.41 ###units are W/m^2




####ok, here we go, here is a MLD description that is going to shift nutrient levels in a seasonal pattern
MLD_ind <- 30*(-cos(seq(0, 2*pi, length = nrow(sun_light)))+1)

plot(sun_time, MLD_ind, type = "o") ###notice that this index has a yearly cycle, just like light (Im manufacturing this using the above cosine function)
###this small loop will give you a sense of the vertical nutrient field is being constructed



z <- -100:0 ######depth vector for our underwater fields. Units are meters


N.field <- NULL #######creating the underwater Nutrient field...the MLD is determined by the index above and follows a -cosine pattern
for(x in 1:length(PAR)){
    
    N.vec <- 0.0248/(1+exp((z+MLD_ind[x])/5))
    
    N.field <- rbind(N.field, (N.vec))
}####end x

image.plot(sun_time, z, N.field) ###notice how the nutrient concentrations for the simulation ##it will be recomputed in the loop below, but this is just to give you an image of it




uI.field <- NULL ###creating the underwater light field (it will be filled with values in the next loop)
pp.field <- NULL ######creating the underwater primary production field (it will be filled with values in the next loop)
N.field <- NULL #######creating the underwater Nutrient field...the MLD is determined by the index above and follows a -cosine pattern
mu.field <- NULL ######this object will hold the interaction relationship between productivity and growth rate
N.mu.field <- NULL ###this object will hold the relationship between N and mu
C.mu.field <- NULL ###this object will hold the relationship between N and mu
cell.field <- NULL #####this object will hold the solution of the growth equation

sum.cell.vec <- NULL ###this object will hold the changing total cell numbers in the water column
k.field <- NULL ######this object will hold the changing mean attenuation coeffiecients as cell numbers change
sum.prod.vec <- NULL ######water column integrated primary production


cell.vec <- c(rep(10000, times = length(z)-1), NA) ####setting the initial cell field that will be changed as cells grow

###take time to understand how the loop is working

for(x in 1:length(PAR)){
    
    N.vec <- 0.0248/(1+exp((z+MLD_ind[x])/5))  ###vertical profile of nutrients mgN/m^3 ###This is turning the MLD_index into a vertical profile of nutrients
    uI.vec <- PAR[x]
    pp.vec <- NA ####the NA is a place holder for the surface
    k.vec <- cell.vec/100000 #####this is a very rough estimation of how cell numbers may impact diffuse attenuation coeffiecient #####attenuation coefficient. Unit is m^-1
    N_Km = 0.06 ##kgN/dy Half saturation constant for computing growth rate (mu)
    N_mumax <- 1 ####units in dy^-1*kgN^-1 Maximum growth rate of phytoplankton
    C_Km = .06 ##mgC/dy Half saturation constant for computing growth rate (mu)
    C_mumax <- 1 ####units in dy^-1 Maximum growth rate of phytoplankton
    DR = -0.1
    
    for(y in 2:length(z)){
        
        pmax = 5 ####units are mgC/dy
        Ek = 50 ####units are W/m2
        R = 1 ###units are mgC/dy
        
        uI <- uI.vec[y-1]*exp(rev(k.vec)[y]*-1)######Lambert-Beer Law solved at every depth (delta z = -1) ##units are W/m2
        uI.vec <- c(uI.vec, uI)#####building the underwater light field with the results of the lambert-Beer Law
        ###NOTE, in previous versions of the mode, we had a stocastic term for each physiological paramter.
        #####Now we have moved the stocastic term into the productivity equation. In this way, we are integrating any uncertainty in the physiology and lumping it into one term that relates to productivity.
        pp <- pmax*tanh(uI/Ek) - R  ###units are mgC/hr ###notice that we have folded the stocatic portion into one equation
        pp.vec <- c(pp.vec, pp)
        
        
    }###end y
    
    uI.vec <- rev(uI.vec)  ###just reversing these vectors because of the way we built them in the previous loop
    pp.vec <- rev(pp.vec)
    
    N_mu <- (N_mumax*N.vec/(N_Km+N.vec)) ########relating growth rate to nitrogen via Michalis-Menton(big assumption here) using a fixed mumax and Km (nobody beleives that they are constant, but this is just a model). We are also adding a loss term (this means that N can be so low that growth rate is negative)
    N_mu <- replace(N_mu, N_mu < 0, DR)
    C_mu <- (C_mumax*abs(pp.vec))/(C_Km+abs(pp.vec)) ########relating growth rate to producivity via Michalis-Menton(big assumption here) using a fixed mumax and Km (nobody beleives that they are constant, but this is just a model).
    #####ok, this is a bit of a short-cut. The this line makes sure that there are negative growth rates where there is negative productivity. Again, just a model for instructive purposes.
    C_mu <- replace(C_mu, pp.vec < 0, DR)
    
    mu.vec <- N_mu*C_mu ####this is how the growth rate determined by N and the the growth rate determined by C are related.
    ######Nobody knows what this is for sure, so it is just a shortcut
    
    mu.vec[which(N_mu<0 | C_mu<0)] <- DR #####this fixes a problem with our above interaction
    
    #####this is also a short cut. This makes sure that growth due to N is 0 if the cell is a respirer.
    
    cell.vec <- cell.vec*exp(mu.vec*1)#####this computes the new population
    
    sum.cell.vec <- c(sum.cell.vec, sum(cell.vec, na.rm = T)) ##this line sums the total number of cells in the water column
    sum.prod.vec <- c(sum.prod.vec, sum(pp.vec, na.rm = T)) ##this line sumes the total amount of productivity in the water column
    
    uI.field <- rbind(uI.field, (uI.vec))
    k.field <- rbind(k.field, (k.vec))
    pp.field <- rbind(pp.field, (pp.vec))
    N.field <- rbind(N.field, (N.vec))
    N.mu.field <- rbind(N.mu.field, N_mu)
    C.mu.field <- rbind(C.mu.field, C_mu)
    mu.field <- rbind(mu.field, (mu.vec))
    cell.field <- rbind(cell.field, (cell.vec))
    print(x)
}#####end x

#####make a color image of the underwater light field
image.plot(sun_time, (z), uI.field, ylab = "Depth", xlab = "Time") ####see ?image.plot

#####make a color image of the underwater primary production field
image.plot(sun_time, (z), pp.field, ylab = "Depth", xlab = "Time") ####see ?image.plot
#####plot the relationship between uI.field and pp.field

plot(uI.field, pp.field, xlab = "Light W/m2", ylab = "Primary Production mgC/hr")



######Question 1: Using the below plot, describe why growth rate for Carbon has the pattern it does. Be through!
image.plot(sun_time, z, C.mu.field)

######Question 2: Using the below plot, describe why growth rate for Nitrogen has the pattern it does.
image.plot(sun_time, z, N.mu.field)

#####Question 3: Describe the pattern of total growth rate in terms of C and N concentrations in the below plot. Why does this pattern emerge?
image.plot(sun_time, z, mu.field)


######Question 4: Describe the pattern of cell number during this simulation. What factors are drivng it?
image.plot(sun_time, z, cell.field)

######Question 5: The below plot is a contour plot of the N field (black) and the cell field (green). Describe the relationship between these contours. What season is cell abundance the hightest?

contour(sun_time, z, N.field, lwd=2)
contour(sun_time, z, cell.field, col = "green", lwd =2, add=T)
