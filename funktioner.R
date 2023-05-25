# Funktioner

# Skapar en uppdelad raster med alla rutor
raster_grid <- function(raster, boxwidth, boxheight, buffer){
  
  #Delar upp raster i delar med storlek boxheight, boxwidth
  
  #R?knar ut h?jd och l?ngd
  height <- ymax(raster) - ymin(raster)
  width <- xmax(raster) - xmin(raster)
  
  #R?knar ut antal rutor (avrundar upp?t)
  hojdled <- ceiling(height/boxheight)
  sidled <- ceiling(width/boxwidth)
  antal_rutor <- hojdled * sidled
  
  #Skapar en 2d-lista att fylla med raster senare
  raster_lista <- as.list(rep(NA, antal_rutor))
  dim(raster_lista) <- c(sidled, hojdled)
  
  #Fyll listan med raster
  message("Delar upp raster...")
  progress <- 1
  start_x <- xmin(raster)
  start_y <- ymin(raster)
  
  for (i in 1:sidled){
    
    for (j in 1:hojdled){
      
      #?ndrar crop beroende p? position av i och j
      
      #Nedre v?nstra h?rnet
      if (i==1 & j==1){
        cropbox <- c(start_x, start_x+(boxwidth+buffer), 
                     start_y, start_y+(boxheight+buffer))
      }
      
      #V?nsterkanten
      else if (i==1){
        cropbox <- c(start_x, start_x+(boxwidth+buffer), 
                     start_y+(j-1)*boxheight-buffer, start_y+(j-1)*boxheight+(boxheight+buffer))
      }
      
      #Nedre kanten
      else if (j==1){
        cropbox <- c(start_x+(i-1)*boxheight-buffer, start_x+(i-1)*boxheight+(boxwidth+buffer), 
                     start_y, start_y+(boxheight+buffer))
      }
      
      #?vre h?gra h?rnet
      else if ( (i+boxwidth+buffer)>xmax(raster) & (j+boxwidth+buffer)>ymax(raster)){
        cropbox <- c(start_x+(i-1)*boxwidth-buffer, xmax(raster), 
                     start_y+(j-1)*boxheight-buffer, ymax(raster))
      }
      
      #H?gerkanten
      else if ( (i+boxwidth+buffer)>xmax(raster)){
        cropbox <- c(start_x+(i-1)*boxwidth-buffer, xmax(raster), 
                     start_y+(j-1)*boxheight-buffer, start_y+(j-1)*boxheight+(boxheight+buffer))
      }
      
      #?vre kanten
      else if ( (j+boxwidth+buffer)>ymax(raster) ){
        cropbox <- c(start_x+(i-1)*boxwidth-buffer, start_x+(i-1)*boxwidth+(boxwidth+buffer), 
                     start_y+(j-1)*boxheight-buffer, ymax(raster))
      }
      
      else{
        cropbox <- c(start_x+(i-1)*boxwidth-buffer, start_x+(i-1)*boxwidth+(boxwidth+buffer), 
                     start_y+(j-1)*boxheight-buffer, start_y+(j-1)*boxheight+(boxheight+buffer))
      }
      
      raster_lista[[i, j]] <- crop(raster, cropbox)
      
      message(progress, "/", antal_rutor)
      progress <- progress + 1
      
    }
    
  }
  print("Raster uppdelad enligt plot.")
  
  #Ritar ut uppdelning som gjordes
  vlinjer <- seq(start_x, xmax(raster), by=boxwidth)
  hlinjer <- seq(start_y, ymax(raster), by=boxheight)
  
  #Bufferlinjer, vertikala
  vbuff_plus <- seq(start_x+boxheight+buffer, xmax(raster), by=boxwidth)
  vbuff_minus <- seq(start_x+(boxheight-buffer), xmax(raster)-boxheight, by=boxheight)
  
  #Bufferlinjer, horisontella
  hbuff_plus <- seq(start_y+boxwidth+buffer, ymax(raster), by=boxwidth)
  hbuff_minus <- seq(start_y+(boxwidth-buffer), ymax(raster)-boxwidth, by=boxwidth)
  
  #Ritar ut
  plot(raster)
  
  if (buffer!=0){
    abline(v=vbuff_plus, col="purple")
    abline(v=vbuff_minus, col="purple")
    abline(h=hbuff_plus, col="purple")
    abline(h=hbuff_minus, col="purple")
  }
  return(raster_lista)
  
}

# Skapar kriging-modeller med individuella variogram med en klustrad raster
krig_grid <- function(raster_grid, vario_typ){
  
  #Skapar en tom lista att fylla med krigade rasters
  antal_celler <- nrow(raster_grid)*ncol(raster_grid)
  kriging_lista <- as.list(rep(NA, antal_celler))
  dim(kriging_lista) <- c(nrow(raster_grid), ncol(raster_grid))
  
  #Skattar variogram f?r varje cell
  variogram_lista <- kriging_lista
  #Skapar en lista att spara spdf f?r varje cell
  spdf <- kriging_lista
  
  message("Tar ut variogram...")
  progress <- 1
  for (i in 1:nrow(variogram_lista)){
    
    for (j in 1:ncol(variogram_lista)){
      
      #G?r om rastercell till spdf
      spdf[[i,j]] <- as.data.frame(raster_grid[[i,j]], xy=TRUE)
      colnames(spdf[[i,j]]) <- c("x", "y", "value")
      coordinates(spdf[[i,j]]) <- ~ x + y
      #Tar bort NA-v?rden
      spdf[[i,j]] <- spdf[[i,j]][is.na(spdf[[i,j]]$value)==FALSE,]
      
      #Anpassar variogram
      lzn.vgm <- variogram(log(value)~1, spdf[[i,j]]) 
      lzn.fit <- fit.variogram(lzn.vgm, model=vgm(vario_typ))
      
      variogram_lista[[i,j]] <- lzn.fit
      message(progress, "/", antal_celler)
      progress <- progress +1
    }
    
  }
  
  #K?r kriging f?r varje cell
  message("Kör kriging...")
  progress <- 1
  for (i in 1:nrow(variogram_lista)){
    
    for (j in 1:ncol(variogram_lista)){
      
      if (nrow(spdf[[i,j]])<length(raster_grid[[i,j]])){
        x.range <- range(spdf[[i,j]]$x)
        y.range <- range(spdf[[i,j]]$y)
        grid <- expand.grid(x = seq(x.range[1], x.range[2], by=2),
                            y = seq(y.range[1], y.range[2], by=2))
        coordinates(grid) <- ~x + y
        
        kriging_lista[[i,j]] <- krige(log(value)~1, spdf[[i,j]], grid, model=variogram_lista[[i,j]])
        kriging_lista[[i,j]] <- as.data.frame(kriging_lista[[i,j]])[,1:3]
        colnames(kriging_lista[[i,j]]) <- c("x", "y", "value")
        
      }else{
        
        kriging_lista[[i,j]] <- as.data.frame(spdf[[i,j]])
        colnames(kriging_lista[[i,j]]) <- c("x", "y", "value")
        kriging_lista[[i,j]]$value <- log(kriging_lista[[i,j]]$value)
        
      }
      
      message(progress, "/", antal_celler)
      progress <- progress +1
    }
  }
  return(list(kriging_lista, variogram_lista))
}

# Skapar kriging-modeller med fixt variogram med en klustrad raster
krig_grid_fixvario <- function(raster_grid, vario_obj){
  
  #Skapar en tom lista att fylla med krigade rasters
  antal_celler <- nrow(raster_grid)*ncol(raster_grid)
  kriging_lista <- as.list(rep(NA, antal_celler))
  dim(kriging_lista) <- c(nrow(raster_grid), ncol(raster_grid))
  
  #Skapar en lista att spara spdf f?r varje cell
  spdf <- kriging_lista
  
  #Kör kriging f?r varje cell
  message("Kör kriging...")
  progress <- 1
  for (i in 1:nrow(spdf)){
    
    for (j in 1:ncol(spdf)){
      #Gör om rastercell till spdf
      spdf[[i,j]] <- as.data.frame(raster_grid[[i,j]], xy=TRUE)
      colnames(spdf[[i,j]]) <- c("x", "y", "value")
      coordinates(spdf[[i,j]]) <- ~ x + y
      #Tar bort NA-v?rden
      spdf[[i,j]] <- spdf[[i,j]][is.na(spdf[[i,j]]$value)==FALSE,]
      
      if (nrow(spdf[[i,j]])<length(raster_grid[[i,j]])){
        x.range <- range(spdf[[i,j]]$x)
        y.range <- range(spdf[[i,j]]$y)
        grid <- expand.grid(x = seq(x.range[1], x.range[2], by=2),
                            y = seq(y.range[1], y.range[2], by=2))
        coordinates(grid) <- ~x + y
        
        kriging_lista[[i,j]] <- krige(log(value)~1, spdf[[i,j]], grid, model=vario_obj)
        kriging_lista[[i,j]] <- as.data.frame(kriging_lista[[i,j]])[,1:3]
        colnames(kriging_lista[[i,j]]) <- c("x", "y", "value")
        
      }else{
        
        kriging_lista[[i,j]] <- as.data.frame(spdf[[i,j]])
        colnames(kriging_lista[[i,j]]) <- c("x", "y", "value")
        kriging_lista[[i,j]]$value <- log(kriging_lista[[i,j]]$value)
        
      }
      
      message(progress, "/", antal_celler)
      progress <- progress +1
    }
  }
  return(list(kriging_lista))
}

# Skapar IDW-modeller med en klustrad raster
IDW_grid <- function(raster_grid){
  
  #Skapar en tom lista att fylla med idw rasters
  antal_celler <- nrow(raster_grid)*ncol(raster_grid)
  IDW_lista <- as.list(rep(NA, antal_celler))
  dim(IDW_lista) <- c(nrow(raster_grid), ncol(raster_grid))
  
  #Skapar en lista att spara spdf för varje cell
  spdf <- IDW_lista
  
  #Kör idw för varje cell
  message("Kör idw...")
  progress <- 1
  for (i in 1:nrow(IDW_lista)){
    
    for (j in 1:ncol(IDW_lista)){
      
      #G?r om rastercell till spdf
      spdf[[i,j]] <- as.data.frame(raster_grid[[i,j]], xy=TRUE)
      colnames(spdf[[i,j]]) <- c("x", "y", "value")
      coordinates(spdf[[i,j]]) <- ~ x + y
      
      #Tar bort NA-v?rden
      spdf[[i,j]] <- spdf[[i,j]][is.na(spdf[[i,j]]$value)==FALSE,]
      
      if (nrow(spdf[[i,j]])<length(raster_grid[[i,j]])){
        x.range <- range(spdf[[i,j]]$x)
        y.range <- range(spdf[[i,j]]$y)
        grid <- expand.grid(x = seq(x.range[1], x.range[2], by=2),
                            y = seq(y.range[1], y.range[2], by=2))
        coordinates(grid) <- ~x + y
        
        IDW_lista[[i,j]] <- idw(log(value)~1, spdf[[i,j]], grid, idp=4)
        IDW_lista[[i,j]] <- as.data.frame(IDW_lista[[i,j]])[,1:3]
        colnames(IDW_lista[[i,j]]) <- c("x", "y", "value")
        
      }else{
        
        IDW_lista[[i,j]] <- as.data.frame(spdf[[i,j]])
        colnames(IDW_lista[[i,j]]) <- c("x", "y", "value")
        IDW_lista[[i,j]]$value <- log(IDW_lista[[i,j]]$value)
        
      }
      
      message(progress, "/", antal_celler)
      progress <- progress +1
    }
  }
  return(list(IDW_lista))
}

spline_grid <- function(raster_grid){
  antal_celler <- nrow(raster_grid)*ncol(raster_grid)
  
  spline_lista <- as.list(rep(NA, antal_celler))
  dim(spline_lista) <- c(nrow(raster_grid), ncol(raster_grid))
  
  df <- spline_lista
  
  message("Kör thin plate spline...")
  progress <- 1
  for (i in 1:nrow(raster_grid)){
    
    for (j in 1:ncol(raster_grid)){
      
      df[[i,j]] <- as.data.frame(raster_grid[[i,j]], xy=TRUE)
      colnames(df[[i,j]]) <- c("x", "y", "value")
      df[[i,j]] <- df[[i,j]][is.na(df[[i,j]]$value)==FALSE,]
      
      if (nrow(df[[i,j]])<length(raster_grid[[i,j]])){
        x.range <- range(df[[i,j]]$x)
        y.range <- range(df[[i,j]]$y)
        grid <- expand.grid(x = seq(x.range[1], x.range[2], by=2),
                            y = seq(y.range[1], y.range[2], by=2))
        grid$coords <- matrix(c(grid$x, grid$y), byrow=F, ncol=2)
        
        df[[i,j]]$coords <- matrix(c(df[[i,j]]$x, df[[i,j]]$y), byrow=F, ncol=2)
        
        #Predikterar med spline
        
        df_surf <-Tps(x = df[[i,j]]$coords, Y = df[[i,j]]$value)
        df_surf_pred <- predict.Krig(df_surf, grid$coords)
        
        spline_lista[[i,j]] <- as.data.frame(cbind(grid[,1:2], "value" = log(df_surf_pred)))
        colnames(spline_lista[[i,j]]) <- c("x", "y", "value")
        
      }else{
        
        spline_lista[[i,j]] <- as.data.frame(df[[i,j]])
        colnames(spline_lista[[i,j]]) <- c("x", "y", "value")
        spline_lista[[i,j]]$value <- log(spline_lista[[i,j]]$value)
        
      }
      
      message(progress, "/", antal_celler)
      progress <- progress +1
    }
  }
  return(list(spline_lista))
  
}

grid_to_raster <- function(krig_grid){
  library(dplyr)
  
  #S?tter ihop alla till en spatial points dataframe
  spdf_grid <- krig_grid[[1]]
  spdf_alla <- data.frame(x=NA, y=NA, value=NA)
  
  for (i in 1:nrow(spdf_grid)){
    for (j in 1:ncol(spdf_grid)){
      spdf_alla <- rbind(spdf_alla, as.data.frame(spdf_grid[i,j]))
    }
  }
  # Tar bort första raden av spdf som inneh?ller NA
  spdf_alla <- spdf_alla[-1,]
  
  # Aggregera över x och y
  spdf_agg <- spdf_alla %>%
    group_by(x, y) %>%
    summarize(value1 = mean(value, na.rm = TRUE))
  spdf_agg <- as.data.frame(spdf_agg)
  
  #Anti-loggar value
  spdf_agg$value1 <- exp(spdf_agg$value1)
  
  raster <- rasterFromXYZ(spdf_agg)
  return(raster)
}



