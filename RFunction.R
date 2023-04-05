library(sf)
library(adehabitatLT)
library(move)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

rFunction <-function(data, barrier_files = NULL, buffer=1000, b_time=4, p_time=36, w=72,
                     max_cross = 0,  sd_multiplier = 1,units = "hours")
{
  Sys.setenv(tz="UTC")
  
  #### Some R functions updated from the BABA project to be used in this 
  #increase movement segment by one points before and one point after the focused encounter ####
  movement.segment.b1 <- function(animal, pt1, pt2) {
    animal <-animal[order(animal$timestamp),]
    animal$ptsID <-seq(nrow(animal))
    
    segments <- animal[animal$ptsID >= pt1 - 1 &  animal$ptsID <= pt2 + 1, ]
    # seg.line <- Lines(Line(coordinates(segments)), ID = segments$timestamp[1])
    # 
    # segments.sp <- SpatialLines(list(seg.line), proj4string = animal@proj4string)
    seg.line <-st_as_sf(segments, coords=c("location.long", "location.lat"), crs=st_crs(4326))
    segments.sf <- seg.line %>% group_by (trackId) %>% 
      summarise(do_union = FALSE) %>%  st_cast("LINESTRING")
    return(segments.sf)
  }
  
  # calculate straigness of movement segment ####
  strtns1 <- function(mov_seg) {
    
    if (sum(duplicated(mov_seg$timestamp)) > 0 ) {
      straightness = NA
      # warning("There are duplicated timestamps")
    } else {
      
      # calculate trajectory
      traj <- adehabitatLT::as.ltraj(xy = st_coordinates(mov_seg), date = mov_seg$timestamp, id = as.character(mov_seg$trackId))
      
      #moving distance from first pt to last pt in the burst
      traj.dist <- sqrt(
        (traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])]) * (traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])]) +
          (traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])]) * (traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])])
      )
      
      #sum of all step lengths
      traj.lgth <- sum(traj[[1]]$dist, na.rm = TRUE)
      
      #straightness ranges from 0 to 1. More close to 0 more sinuous it is.
      straightness <- traj.dist/traj.lgth }
    
    return(straightness)
  }
  
  ### Loading the road data and cropping according to the bounding box of the data
  #roads <- st_read("./data/GRIP_roads_NASAY2Y/GRIP_roads_NASAY2Y.shp")
  roads <- st_read(paste0(getAppFilePath("barrier_files"),"roads.shp"))
  roads_crop <- st_crop(roads, st_bbox(data))
  roads_buffer <-st_buffer(roads_crop, dist= buffer)
  
  data_df <- as.data.frame(data)
  coo <- data.frame(coordinates(data))
  names(coo) <- c("location.long","location.lat")
  names(data_df) <- make.names(names(data_df),allow_=FALSE)
  data_df <- data.frame(data_df,coo)
  data_df$timestamp<- move::timestamps(data)
  data_df$trackId<- move::trackId(data)
  data_df$ptsID <-NA
  
  uid <-unique(data_df$trackId)
  interval_dat <- data.frame(uid)
  rd_enctr <-list()
  
  for(i in 1:length(uid))
  {
    temp_dat <- subset(data_df, data_df$trackId == uid[i])
    temp_dat <-temp_dat[order(temp_dat$timestamp),]
    temp_dat$ptsID <-seq(nrow(temp_dat))
    
    interval_dat$fix_int[i] <- median(as.numeric(diff(temp_dat$timestamp), units = "hours"))
    class(interval_dat$fix_int)<-"numeric"
    
    temp_dat_sf <-st_as_sf(temp_dat,coords=c("location.long", "location.lat"), 
                           crs= crs(roads_buffer))
    rd_enctr_temp <-st_intersection(temp_dat_sf, roads_buffer)
    rd_enctr_temp <-rd_enctr_temp[,names(rd_enctr_temp) %in% names(temp_dat_sf)]
    
    if (nrow(rd_enctr_temp) == 0) {
      warning(paste0 ("Individual ", uid[i], " has no locations overlapped with the barrier buffer and is eliminated from analysis." ))
      next()
    }
    
    rd_enctr_temp$timediff <- c(interval_dat$fix_int[i], as.numeric(diff(rd_enctr_temp$timestamp), units = "hours"))
    rd_enctr_temp$ptdiff <- c(2, as.numeric(diff(rd_enctr_temp$ptsID)))
    
    rd_enctr_temp$timediff2 <- round(rd_enctr_temp$timediff - interval_dat$fix_int[i], digits = 1)
     
    ## then, if any timediff2 is > interval but <= tolerance, we need to bring in the missing points from outside the buffer.
    if(any(rd_enctr_temp$ptdiff > 1)) {
      
      idx_pts_of_interest <- which(rd_enctr_temp$ptdiff > 1)
      
      for(pt in idx_pts_of_interest) {
        # find out what pts to fetch
        ptsID_of_interest_B <- rd_enctr_temp$ptsID[pt]+1
        ptsID_of_interest_A <- rd_enctr_temp$ptsID[pt]-1
        
        # fetch the points outside of the buffer and placehold timediff as NA and timediff2 as 0
        fetched_pt <- rbind(temp_dat_sf[ temp_dat_sf$ptsID == ptsID_of_interest_A,] , 
                            temp_dat_sf[ temp_dat_sf$ptsID == ptsID_of_interest_B,])
        
        if (nrow(fetched_pt) == 0) {  # if there's no point outside of the buffer between the timestamp that means there's missing data
          # since the missing data is still within the tolerence, we consider timediff2=0 so the points before and after will be in the same event
          rd_enctr_temp$timediff2[pt] <- 0
          next() } 
        else {
          fetched_pt$timediff <- NA
          fetched_pt$timediff2 <- 0 
          fetched_pt$ptdiff <- rd_enctr_temp$ptdiff[pt]
          # replace timediff2 of pts_of_interests to 0
          rd_enctr_temp$timediff2[pt] <- 0 
          # append fetched points to each other 
          if(pt == idx_pts_of_interest[1]) fetched_pts <- fetched_pt 
          else fetched_pts <- rbind(fetched_pts, fetched_pt)
        }
      }
      
      # append fetched pts
      rd_enctr_temp <- rbind(rd_enctr_temp, fetched_pts)
      #recorder animal i encounter event dataframe
      rd_enctr_temp <- rd_enctr_temp[order(rd_enctr_temp$ptsID), ]
    }
    
    ## then do the cum sum of the new dataframe based on timediff2, using that as the unique burst ID (with animalID) 
    rd_enctr_temp$burstID <- paste(uid[i], rd_enctr_temp$ptdiff, sep = "_")
    
    interval_dat[i,3] <- nrow(rd_enctr_temp)
    rd_enctr[[i]] <-rd_enctr_temp
    
    }
  
    encounter_dat <-do.call(rbind, rd_enctr)  
  
    rs <-raster(extent(data), res = .05, crs= projection(roads)) ## create the raster of the dataset
    enc_dat <-cbind(encounter_dat$coords.x1, encounter_dat$coords.x2) ##extract the coordinates
    
    tab <- table(cellFromXY(rs, enc_dat))
    rs[as.numeric(names(tab))] <- tab
    d <- data.frame(coordinates(rs), count=rs[])
    d<-d[complete.cases(d$count),]
    
  ### Plot the encounter locations of road and the buffer  
density_plot <-  ggplot()+geom_sf(data=roads_buffer, size=0.5)+
    geom_sf(data=roads_crop, col="brown", size=1)+
    geom_tile(data= d, aes(x= x, y=y, fill = count), alpha=0.9)+
    scale_fill_distiller(palette = "Spectral", direction = -1)+
    labs(x= "Longitude", y= "Latitude", fill= "Number of animal\nlocations")+
    #geom_sf(data=encounter_dat) + 
    theme_bw()

  ##########################################################################################
  ############ Classifying barrier behaviors
  
  ### create empty object that will hold results ----
  event_df <- NULL
  plot_list <-list()
  
  ## Classifying short events ##
  ## run classification procedure for each encounter ####
  for(i in unique(encounter_dat$burstID)) {
    
    # get what we need from the encounter ####
    encounter_i <- encounter_dat[encounter_dat$burstID == i, ]
    animal_i <- data_df[data_df$trackId == encounter_i$trackId[1],]
    start_time <- encounter_i$timestamp[1]
    end_time <- encounter_i$timestamp[nrow(encounter_i)]
    duration <-  difftime (end_time, start_time, units = "hours")
    
    # calculating straightness of the encounter event ###
    ## this will be used for median duration event but is output for reference for other event ####
    straightness_i <- strtns1(encounter_i)
    
    # classify short encounters (bounce and quick cross) ####
    # if no more than b*interval, only spend small amount of time in this burst
    if (duration <= b_time) {
      pt.first <- encounter_i$ptsID[1]#first point in the burst
      pt.last <- encounter_i$ptsID[nrow(encounter_i)]
      
      # extract movement segment with one point before and one point after the segmentation ####
      mov_seg_i <- movement.segment.b1(animal_i, pt.first, pt.last)
      
      # count the number of crossing ####
      int.num <- nrow(st_intersection(mov_seg_i, roads_crop))
      
      # if no crossing and we didn't have both points (before and after), then we can't tell if it crossed
      if (int.num == 0 & nrow(st_coordinates(mov_seg_i)) != (nrow(encounter_i)+2)) {
        # means that no points were before or after the encounter and we can't tell if the animal crossed
        classification <- "unknown"
        
      } else {
        classification <- ifelse(int.num == 0, "Bounce", "Quick_Cross")
      }
      
      # plot these examples to check later
      
      plot_list[[i]]<- ggplot()+
        geom_sf(data = roads_buffer, size=0.5)+
        geom_sf(data = roads_crop, size=1, col= "brown")+
        geom_sf(data = mov_seg_i)+
        geom_sf(data = encounter_i, size=1, col ="blue")+
        lims(x=c(st_bbox(mov_seg_i)[1]-0.05,st_bbox(mov_seg_i)[3]+0.05),
             y=c(st_bbox(mov_seg_i)[2]-0.05,st_bbox(mov_seg_i)[4]+0.05))+
        labs(title = paste(i, "_",classification))+
        theme_bw()
        
    }
  
    if (duration > b_time) {
      
      ## first calculate number of crossings (without looking at extra points like we did for short encounter)
      #mov_seg_i <- SpatialLines(list(Lines(Line(st_coordinates(encounter_i)), ID = encounter_i$date[1])), proj4string = animal@proj4string)
      mov_seg_i <- encounter_i %>% group_by (trackId) %>% 
        summarise(do_union = FALSE) %>% st_cast("LINESTRING")
      int.num <- nrow(st_intersection(mov_seg_i, roads_crop))
      ## then check if duration is smaller of bigger than p and classify accordingly
      if(duration > p_time) {
        
        classification <- "Trapped"
        
      } else {
        
        classification <- "TBD" # these will be further classified in the next loop
        
      }
      
      
      # plot these examples to check later
      if (!classification %in% "TBD") {
        plot_list[[i]]<- ggplot()+
          geom_sf(data = roads_buffer, size=0.5)+
          geom_sf(data = roads_crop, size=1, col= "brown")+
          geom_sf(data = mov_seg_i)+
          geom_sf(data = encounter_i, size=1, col ="blue")+
          lims(x=c(st_bbox(mov_seg_i)[1]-0.05,st_bbox(mov_seg_i)[3]+0.05),
               y=c(st_bbox(mov_seg_i)[2]-0.05,st_bbox(mov_seg_i)[4]+0.05))+
          labs(title = paste(i, "_",classification))+
          theme_bw()
      }
    }
  
    event_df <- rbind(event_df, data.frame(Individual_ID = encounter_i$trackId[1],
                                           burstID = i,
                                           long = ifelse(nrow(encounter_i)>1 , st_coordinates(encounter_i)[2,1],st_coordinates(encounter_i)[1,1] ),
                                           lat  = ifelse(nrow(encounter_i)>1 , st_coordinates(encounter_i)[2,2],st_coordinates(encounter_i)[1,2] ),
                                           start_time,    end_time,   duration,
                                           cross = int.num,
                                           straightness =  straightness_i,
                                           eventTYPE = classification,
                                           stringsAsFactors = F))
    }
  
   #############################################
  ###### classify Back-n-forth and Trace ######
  #############################################
  
  ### back-n-forth and trace are based on comparing average straightness around the encounter event
  
  for (i in 1:nrow(event_df)) {
    if (event_df[i, ]$eventTYPE == "TBD") {
      event_i <- event_df[i, ]
      duration_i <- event_i$duration
      straightness_i <- event_i$straightness
      interval <-interval_dat$fix_int[uid == event_i$Individual_ID]
      
      # get movement data of the individual of interest
      animal_i <- data_df[data_df$trackId == event_i$Individual_ID, ]
      encounter_i <- encounter_dat[encounter_dat$burstID %in% event_i$burstID,]
      
      # else {
      #   animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Animal.ID == event_i$AnimalID] & encounter$burstID %in% event_i$burstID, ]
      # }
      
      # keep only data X units before and X after event
      animal_i <- animal_i[animal_i$timestamp >= event_i$start_time - as.difftime(w/2, units = "hours") & 
                           animal_i$timestamp <= event_i$end_time +  as.difftime(w/2, units = "hours"), ]
      
      # identify continuous sections in the remaining movement data
      #animal_i$continuousID <- cumsum(abs(c(interval, round(as.numeric(diff(animal_i$timestamp), units = "hours"), digits = 1)) - interval)) # sep 11, 2020. added abs() to accomodate potential data points with smaller time intervals
      animal_i$continuousID <- animal_i$trackId
        
      # for each continuous sections, calculate straightness of all movements lasting the duration of our event (moving window of the size of the encounter)
      straightnesses_i <- NULL
      for(ii in unique(animal_i$continuousID)) {
        animal_ii <- animal_i[animal_i$continuousID == ii, ]
        
        #duration of period
        duration_ii <- difftime(animal_ii$timestamp[nrow(animal_ii)], animal_ii$timestamp[1], units = "hours")
        
        # calculate straightness only if at least as long as encounter event
        if(duration_ii >= duration_i) {
          for(iii in c(1: (which(animal_ii$timestamp > (animal_ii$timestamp[nrow(animal_ii)] - as.difftime(duration_i, units = "hours")))[1]-1))) {
            mov_seg <- animal_ii[iii: min((iii + as.numeric(duration_i/interval)), nrow(animal_ii)), ]
            mov_seg_sf <-st_as_sf(mov_seg, coords = c("location.long", "location.lat"), crs=st_crs(4326))
            straightnesses_i <- c(straightnesses_i, strtns1(mov_seg_sf))
          }
        }
      }
      
      # make sure there are enough data to calculate average monthly straightness before and after the encounter event
      # (w/interval + 1) is the total possible segments if all data are present. 
      # We define "enough" as at least 1/4 of the total possible segments are present to calculate average straightness.
      if (length(straightnesses_i) >= (w/interval + 1)/4) {
        # minimum max number possible/2 to calculate sd
        upper <- ifelse(is.na(sd(straightnesses_i, na.rm=T)) == FALSE, mean(straightnesses_i, na.rm=T) + sd_multiplier * sd(straightnesses_i, na.rm=T),
                        mean(straightnesses_i) + sd_multiplier *0.01)
        lower <- ifelse(is.na(sd(straightnesses_i, na.rm=T)) == FALSE, mean(straightnesses_i, na.rm=T) - sd_multiplier * sd(straightnesses_i, na.rm=T),
                        mean(straightnesses_i) - sd_multiplier *0.01)
        
        if(straightness_i < lower) event_df[i, ]$eventTYPE <- ifelse(event_i$cross < max_cross, "Back_n_forth", "unknown")
        if(straightness_i > upper) event_df[i, ]$eventTYPE <- ifelse(event_i$cross < max_cross, "Trace", "unknown")
        if(straightness_i >= lower & event_i$straightness <= upper) event_df[i, ]$eventTYPE <- "Average_Movement"
        
      } else {
        event_df[i, ]$eventTYPE = "unknown"
        if(is.null(straightnesses_i)) {straightnesses_i <- NA} # adding this to avoid warning message when ploting.
      }
      
      if (classification %in% c("Back_n_forth", "Trace", "Average_Movement")) {
                 plot_list[[i]]<- ggplot()+
                     geom_sf(data = roads_buffer, size=0.5)+
                     geom_sf(data = roads_crop, size=1, col= "brown")+
                     geom_sf(data = mov_seg_i)+
                     geom_sf(data = encounter_i, size=1, col ="blue")+
                     lims(x=c(st_bbox(mov_seg_i)[1]-0.05,st_bbox(mov_seg_i)[3]+0.05),
                          y=c(st_bbox(mov_seg_i)[2]-0.05,st_bbox(mov_seg_i)[4]+0.05))+
                     labs(title = paste(i, "_",classification))+ theme_bw()
      }
    }
  }
  
  ## clean the encounter spdataframe ##
  encounter_data <- encounter_dat[!duplicated(encounter_dat$burstID),]
  encounter_data <- encounter_data[,c("trackId","burstID","timestamp")]
  encounter_data <- merge(encounter_data, event_df[,c("burstID","eventTYPE")])
  
  ### Print all the classified encounters in one pdf file
  pdf(paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Event_plot_output.pdf"))
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  for (i in unique(encounter_dat$burstID)) {
    print(plot_list[[i]])
  }
  dev.off()
  
  ggsave(density_plot, filename = paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Point_count_density.jpeg"),
         width = 9, height = 6, dpi=300, units = "in")
  # #write.csv(encounter_data, file= paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Encounter_data.csv"))
  write.csv(event_df, file= paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Encounter_event_data.csv"))
  
  return(data)
  
}

  
  
 