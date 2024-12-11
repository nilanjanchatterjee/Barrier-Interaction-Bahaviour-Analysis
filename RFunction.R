library(sf)
library(adehabitatLT)
library(move2)
library(ggplot2)
library(dplyr)
library(raster)
library(RColorBrewer)

rFunction <-  function(data, barrier_files = NULL, buffer=1000,  b_time = 4, p_time = 36, w = 72,
           tolerance = 0, units = "hours", max_cross = 0,  sd_multiplier = 1,den_ras = 500,
           round_fixes = F,exclude_buffer =FALSE) 
  {
    
    # initial checks ----------------------------------------------------------
   
    # ## prepare parameters and check input
    # if (class(animal)[1] != "sf") stop("animal needs to be an sf object")
    # if (sf::st_geometry_type(animal)[1] != 'POINT') stop("animal needs to have a POINT geometry")
    # if (class(barrier)[1] != "sf") stop("barrier needs to be an sf object")
    # if (!(sf::st_geometry_type(barrier)[1] %in% c('MULTILINESTRING', 'LINESTRING'))) stop("barrier needs to have either a LINESTRING or MULTILINESTRING geometry")
    # if (!"date" %in% names(animal)) stop("please rename the date column to 'date'")
    # if (!"Animal.ID" %in% names(animal)) stop("please rename the individual ID column to 'Animal.ID'")
    # if (!(inherits(animal$date, "POSIXct"))) stop("date needs to be 'POSIXct' format")
    # if (sum(is.na(animal$date)) > 0) stop("please exclude rows where date is NA")
    # 
    # if(round_fixes){
    #   interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(round(as.numeric(diff(x), units = units),0)))))
    # } else {
    #   interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(as.numeric(diff(x), units = units)))))
    # }
    # if(is.null(interval)) { ## figure out interval (as the most frequent difference in timestamp) if not provided but give an error if not the same for all individuals
    #   if(all(interval_per_individual == interval_per_individual[1])) interval <- as.numeric(interval_per_individual[1]) else stop("Not all individuals have been sampled at the same frequency. Run individuals with different intervals seperately, or double-check whether your date column is cleaned.")
    # } else {
    #   if (any(as.numeric(interval_per_individual) > interval, na.rm = T)) stop("BaBA interval needs to be no smaller than the actual data interval. Also double-check whether your date column is cleaned.") 
    # }
    
    # b <- b_time / interval
    # if(b < 1) stop("interval needs to be set no bigger than b_time")
    # if (round(b) != b) stop("b_time must be divisible by interval")
    # p <- p_time / interval
    # if (round(p) != p) stop("p_time must be divisible by interval")
    # 
    
    roads <- st_read(paste0(getAuxiliaryFilePath("barrier_files"),"roads.shp"))
    roads_crop <- st_crop(roads, st_bbox(data))
    #roads_buffer <-st_buffer(roads_crop, dist= buffer)
    
    
    data <- data |> mutate(location.long = sf::st_coordinates(data)[,1],
                           location.lat = sf::st_coordinates(data)[,2],
                           trackId = mt_track_id(data),
                           timestamp = mt_time(data))
    
    # classification step 1: generate encounter event data.frame -------------------------------------------------
    
    ## create point ID by individual
    data <- data %>% 
      dplyr::arrange(trackId, timestamp) %>% 
      dplyr::group_by(trackId) %>% 
      dplyr::mutate(ptsID = 1:dplyr::n()) %>% 
      dplyr::ungroup()
    
    ## explicitly suppress constant geometry assumption warning by confirming attribute is constant throughout the geometry. See https://github.com/r-spatial/sf/issues/406 for details.
    sf::st_agr(data) <- 'constant'   
    sf::st_agr(roads_crop) <- 'constant'
    
    ## create buffer around barrier
    print("locating encounter events...")
    barrier_buffer <- roads_crop %>% 
      sf::st_buffer(dist = buffer, nQuadSegs = 5) #%>%   ## Note that nQuadSegs is set to 5 as this was the default value for rgeos::gBuffer in previous versions of BaBA
      #sf::st_union()
    
    ## extract points that fall inside the buffer
    encounter <- sf::st_intersection(data, barrier_buffer)
    
    if (nrow(encounter) == 0) stop("no barrier encounter detected.")
    
    ## create unique burstIDs
    for(i in unique(encounter$trackId)){
      if (nrow(encounter %>% dplyr::filter(trackId == i)) == 0) {
        warning(paste0 ("Individual ", i, " has no locations overlapped with the barrier buffer and is eliminated from analysis." ))
        next()
      }
      encounter_i <-
        encounter %>% 
        dplyr::filter(trackId == i) 
        
        temp_interval <- median(as.numeric(diff(encounter_i$timestamp), units = "hours"))
      class(temp_interval)<-"numeric"
      
        ## add time difference
      encounter_i <-encounter_i %>% dplyr::mutate(
          ## time difference between all points in the buffer
          timediff = c(temp_interval, as.numeric(diff(timestamp), units = units)),
          ## remove the interval from that so when there is no missing point, timediff2 should be 0. If <0, replicated timestamp; if >0, missing timestamp
          timediff2 = round(timediff - temp_interval, digits = 1))
      
      ## if any timediff2 is >= interval but <= tolerance, bring in the missing points from outside the buffer
      if(any(encounter_i$timediff2 >= temp_interval & encounter_i$timediff2 <= tolerance, na.rm = T )) {
        idx_pts_of_interest <- which(encounter_i$timediff2 >= temp_interval & encounter_i$timediff2 <= tolerance)
        
        for(pt in idx_pts_of_interest) {
          ## find out what pts to fetch
          ptsID_of_interest_B <- encounter_i$ptsID[pt]
          ptsID_of_interest_A <- encounter_i$ptsID[pt-1]
          
          ## fetch the points outside of the buffer and placehold timediff as NA and timediff2 as 0
          fetched_pt <- 
            data %>% 
            dplyr::filter(trackId == i & 
                            ptsID > ptsID_of_interest_A & 
                            ptsID < ptsID_of_interest_B)
          
          if (nrow(fetched_pt) == 0) {  ## if there's no point outside of the buffer between the timestamp that means there's missing data
            ## since the missing data is still within the tolerance, we consider timediff2=0 so the points before and after will be in the same event
            encounter_i$timediff2[pt] <- 0
            next() } 
          else {
            fetched_pt$timediff <- NA
            fetched_pt$timediff2 <- 0 
            fetched_pt$burstID <- encounter_i[encounter_i$ptsID == ptsID_of_interest_A,]$burstID
            ## reset timediff2 of pts_of_interests to 0
            encounter_i$timediff2[pt] <- 0 
            ## append fetched points to each other 
            if(pt == idx_pts_of_interest[1]) {fetched_pts <- fetched_pt} else if (exists("fetched_pts")) { fetched_pts <- rbind(fetched_pts, fetched_pt) } else {fetched_pts <- fetched_pt}
          }
        }
        
        ## append fetched pts
        encounter_i <- rbind(encounter_i, fetched_pts)
        ## recorder animal i's encounter event data.frame
        encounter_i <- encounter_i[order(encounter_i$ptsID), ]
      }
      
      ## do the cumulative sum of the new data.frame based on timediff2, using that as the updated unique burst ID (with animalID) 
      encounter_i$burstID <- paste(i, cumsum(encounter_i$timediff2), sep = "_")
      
      ## save into encounter_complete
      if(i == unique(encounter$trackId[1])) encounter_complete <- encounter_i else encounter_complete <- rbind(encounter_complete, encounter_i)
    }
    
    encounter <- encounter_complete ## save back as encounter (encounter_complete is bigger as it includes extra points that are within tolerance)
    
    ###create the density plot for the road encounters  
    rs <-raster(extent(as.vector(st_bbox(data))[c(1, 3, 2, 4)]), 
              res = den_ras/100000, crs= projection(roads)) ## create the raster of the dataset
    enc_dat <-cbind(encounter$location.long, encounter$location.lat) ##extract the coordinates
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
    # classification step 2: classify events -------------------------------------------------
    
    print("classifying behaviors...") 

    ## create empty object that will hold results
    event_df <- NULL
    plot_list <- list()
    
    ## run classification procedure for each encounter
    for(i in unique(encounter$burstID)) {
       
      ## get what we need from the encounter
      encounter_i <- encounter[encounter$burstID == i, ]
      animal_i <- data[data$trackId == encounter_i$trackId[1],]
      start_time <- encounter_i$timestamp[1]
      end_time <- encounter_i$timestamp[nrow(encounter_i)]
      duration <-  difftime (end_time, start_time, units = units)
      
      ## calculating straightness of the encounter event
      ## this will be used for median duration events but is output for reference for other events
      straightness_i <- strtns(encounter_i)
      
      #print(nrow(encounter_i))
      ### classify short events (bounce and quick cross) ---------------------------------------------------
      
      ## if no more than b*interval, only spend small amount of time in this burst
      if (duration <= b_time) {
        pt.first <- encounter_i$ptsID[1] ## first point in the burst
        pt.last <- encounter_i$ptsID[nrow(encounter_i)]
        
        ## extract movement segment with one point before and one point after the segmentation
        mov_seg_i <- movement.segment.b(animal_i, pt.first, pt.last)
        
        ## count the number of crossings
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(roads_crop) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## if no crossing is indicated and both before and after points were missing then we cannot tell if the animal crossed
        if (int.num == 0 & nrow(sf::st_coordinates(mov_seg_i)) != (nrow(encounter_i)+2)) {
          classification <- "unknown"
        } else {
          ## if there was not a crossing, classify as bounce, otherwise quick cross
          classification <- ifelse(int.num == 0, "Bounce", "Quick_Cross")
        }
        plot_list[[i]]<- ggplot()+
          geom_sf(data = barrier_buffer, size=0.5)+
          geom_sf(data = roads_crop, size=1, col= "brown")+
          geom_sf(data = mov_seg_i)+
          geom_sf(data = encounter_i, size=1, col ="dodgerblue")+
          geom_sf(data = encounter_i[1,], size=1.5, col ="forestgreen")+
          #geom_sf(data = encounter_i[nrow(encounter),], size=1.5, col ="tomato")+
          lims(x=c(st_bbox(mov_seg_i)[1]-0.05,st_bbox(mov_seg_i)[3]+0.05),
               y=c(st_bbox(mov_seg_i)[2]-0.05,st_bbox(mov_seg_i)[4]+0.05))+
          labs(title = paste(i, "_",classification),
               subtitle = paste("start_time=", encounter_i$timestamp[1], "end_time=", encounter_i$timestamp[nrow(encounter_i)]))+
          theme_bw()
      }
      
      ## dummy variable to ensure desired plotting and output
      tbd.plot <- 0
      
      if (duration > b_time) {
        
        ### classify trapped events -------------------------------------------------
        
        ## first calculate number of crossings (without looking at extra points like we did for short encounter)
        mov_seg_i <- 
          encounter_i %>% 
          dplyr::summarize(do_union = FALSE) %>% 
          sf::st_cast(to = 'LINESTRING')
        
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(roads_crop) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## check if duration is smaller of bigger than p and classify accordingly
        if(duration > p_time) {
          classification <- "Trapped"
        } else {
          classification <- "TBD" ## these will be further classified in the next loop
        }
        
        
        ### classify back-n-forth, trace, and average movement -----------------------------------------
        
        ## process the "TBD" event types as back-n-forth, trace, or average movement
        ## back-n-forth and trace are based on comparing average straightness around the encounter event
        if(classification == 'TBD'){
          
          tbd.plot <- 1
          
          ## remove points that are inside the buffer if user said so
          if (exclude_buffer) {
            animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$trackId == animal_i$trackId[1]], ]
          }
          
          ## keep only data w/2 units before and w/2 after event
          animal_i <- animal_i[animal_i$timestamp >= start_time - as.difftime(w/2, units = units) & animal_i$timestamp <= end_time +  as.difftime(w/2, units = units), ]
          
          ## identify continuous sections in the remaining movement data
          animal_i$continuousID <- cumsum(abs(c(interval, round(as.numeric(difftime(animal_i$timestamp, 
                                                                          lag(animal_i$timestamp, 1)), 
                                                                 units = units), digits = 1)[-1] - interval))) # sep 11, 2020. added abs() to accomodate potential data points with smaller time intervals
          
          ## for each continuous sections, calculate straightness of all movements lasting the duration of our event (moving window of the size of the encounter)
          straightnesses_i <- NULL
          for(ii in unique(animal_i$continuousID)) {
            animal_ii <- animal_i[animal_i$continuousID == ii, ]
            
            ## duration of period
            duration_ii <- difftime(animal_ii$timestamp[nrow(animal_ii)], animal_ii$timestamp[1], units = units)
            
            ## calculate straightness only if at least as long as encounter event
            if(duration_ii >= duration) {
              for(iii in 1:(which(animal_ii$timestamp > (animal_ii$timestamp[nrow(animal_ii)] - as.difftime(duration, units = units)))[1] -1)) {
                mov_seg <- animal_ii[iii:(iii + duration/interval), ]
                straightnesses_i <- c(straightnesses_i, strtns(mov_seg))
              }
            }
          }
          
          ## make sure there are enough data to calculate average straightness before and after the encounter event
          ## (w/interval + 1) is the total possible segments if all data are present. 
          ## We define "enough" as at least 1/4 of the total possible segments are present to calculate average straightness.
          if (length(straightnesses_i) >= (w/interval + 1)/4) {
            ## minimum max number possible/2 to calculate sd
            upper <- mean(straightnesses_i) + sd_multiplier * stats::sd(straightnesses_i)
            lower <- mean(straightnesses_i) - sd_multiplier * stats::sd(straightnesses_i)
            if(straightness_i < lower) classification <- ifelse(int.num <= max_cross, "Back_n_forth", "unknown")
            if (straightness_i > upper) classification <- ifelse(int.num <= max_cross, "Trace", "unknown")
            if(straightness_i >= lower & straightness_i <= upper) classification <- "Average_Movement"
          } else {
            classification <- "unknown"
            if(is.null(straightnesses_i)) {straightnesses_i <- NA} ## add this to avoid warning message when plotting
          }
          plot_list[[i]]<- ggplot()+
            geom_sf(data = barrier_buffer, size=0.5)+
            geom_sf(data = roads_crop, size=1, col= "brown")+
            geom_sf(data = mov_seg_i)+
            geom_sf(data = encounter_i, size=1, col ="dodgerblue")+
            geom_sf(data = encounter_i[1,], size=1.5, col ="forestgreen")+
            #geom_sf(data = encounter_i[nrow(encounter),], size=1.5, col ="tomato")+
            lims(x=c(st_bbox(mov_seg_i)[1]-500,st_bbox(mov_seg_i)[3]+500),
                 y=c(st_bbox(mov_seg_i)[2]-500,st_bbox(mov_seg_i)[4]+500))+
            labs(title = paste(i, "_",classification),
                 subtitle = paste("start_time=", encounter_i$timestamp[1], "end_time=", encounter_i$timestamp[nrow(encounter_i)]))+ 
            theme_bw()
        }
      }
      
      
      ### Consolidate outputs -----------------------------------------------------
      
      ## plot the encounters to check later, if desired

      
      ## combine output
      event_df <- rbind(event_df, data.frame(
        AnimalID = encounter_i$trackId[1],
        burstID = i,
        easting = sf::st_coordinates(encounter_i)[1, 1],
        northing = sf::st_coordinates(encounter_i)[1, 2],
        start_time,
        end_time,
        duration,
        cross = int.num,
        str_i = straightness_i,
        str_mean = ifelse(tbd.plot == 0, NA, mean(straightnesses_i)),
        str_sd = ifelse(tbd.plot == 0, NA, stats::sd(straightnesses_i)),
        eventTYPE = classification,
        stringsAsFactors = F
      ))
    }
    
    
    # finalize data -----------------------------------------------------------
    
    print("creating dataframe...")
    ## clean the encounter data
    encounter_final <- encounter %>% 
      dplyr::filter(!duplicated(burstID)) %>% 
      dplyr::left_join(event_df %>% 
                         dplyr::select(burstID, eventTYPE),
                       by = 'burstID') %>% 
      dplyr::select(trackId, burstID, timestamp, eventTYPE)
    
    ### Print all the classified encounters in one pdf file
    pdf(appArtifactPath("Event_plot_output.pdf"))
    par(mfrow=c(2,2), mar=c(4,4,3,1))
    
    #encounter_final_sub<- subset(encounter_final ,!(eventTYPE %in% c("Bounce", "Quick_Cross")))
    for (i in unique(encounter_final$burstID)) {
      print(plot_list[[i]])
    }
    dev.off()
    
    ## Saving the outputs
    ggsave(density_plot, filename = appArtifactPath("Point_count_density.jpeg"),
           width = 9, height = 6, dpi=300, units = "in")
    # #write.csv(encounter_data, file= paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Encounter_data.csv"))
    write.csv(event_df, file= appArtifactPath("Encounter_event_data.csv"))
    
    return(data)
  
  }


## increase movement segment by one points before and one point after the focused encounter
movement.segment.b <- function(animal, pt1, pt2) {
  pts_tmp <- animal[animal$ptsID >= pt1 - 1 & animal$ptsID <= pt2 + 1, ]
  pts_comb <- dplyr::summarize(pts_tmp, do_union = FALSE)
  segments_out <- sf::st_cast(pts_comb, to = 'LINESTRING')
  return(segments_out)
}

## helper function on for calculating Euclidean distance
calc_dist <- function(x.start, x.end, y.start, y.end){
  sqrt((x.end - x.start)^2 + (y.end - y.start)^2)
}

## calculate straightness of movement segment
strtns <- function(mov_seg) {
  
  locs_tmp <- sf::st_coordinates(mov_seg)
  
  if (sum(duplicated(mov_seg$timestamp)) > 0 ) {
    straightness <- NA
    # warning("There are duplicated timestamps")
  } else if(nrow(locs_tmp) == 1){
    straightness <- NA
  } else {
    ## calculate Euclidean distance between the first and last point
    euc_dist <- 
      as.numeric(
        calc_dist(x.start = locs_tmp[1,1], x.end = locs_tmp[nrow(locs_tmp),1],
                  y.start = locs_tmp[1,2], y.end = locs_tmp[nrow(locs_tmp),2]))
    
    ## calculate path distance as the sum of all step lengths
    mov_seg$dist <- NA
    for(j in 2:nrow(mov_seg)){
      mov_seg$dist[j] <- calc_dist(x.start = locs_tmp[j - 1, 1],
                                   x.end = locs_tmp[j, 1],
                                   y.start = locs_tmp[j - 1, 2],
                                   y.end = locs_tmp[j, 2])
    }
    path_dist <- sum(mov_seg$dist, na.rm = TRUE)
    
    ## calculate straightness as the ratio of Euclidean to path distance.
    ## straightness ranges from 0 to 1 with values closer to 0 being more
    ## sinuous.
    straightness <- euc_dist/path_dist
  }
  
  return(straightness)
}
