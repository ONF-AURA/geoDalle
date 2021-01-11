

#' Raster du potentiel de réserve profondes
#'
#'
#'
#' @param mnt0 raster MNT, à résolution minimale de 5m
#' @param zn_prf polygone (sf) de la zone à étudier
#' @param zn polygone sf de découpage du MNT
#' @param geo_meuble polygones des roches meubles (moraines, molasses, éboulis, alluvions...)
#' @param p.pendage points d'observation BRGM des pendages
#' @param con_geo connexion database
#' @param db_geo nom de la table géologie
#' @param db_pendage nom de la table des points de pendafe
#' @param geo_not_dalle vecteur des premières lettres des notations des couches de roches non meubles
#' @param buff buffer (en m) à utiliser autour de zn_prf pour réaliser les calculs (défaut: 500m)
#'
#' @return raster du risque de présence d'une dalle compact
#'
#' @import magrittr
#' @importFrom raster raster values terrain crop aggregate stack focalWeight focal mask clump extract xFromCell yFromCell
#' @importFrom sf st_transform st_crs read_sf st_union st_buffer as_Spatial
#' @importFrom purrr map
#' @importFrom dplyr group_by select filter pull summarise count mutate arrange desc
#' @importFrom tmap tm_shape tm_text
#' @importFrom spatialEco curvature tri classBreaks
#' @export
#'
#' @example
#' dir <- system.file(package = "geoDalle")
#' mnt <- raster(file.path(dir, "data","mnt0.tif"))
#' zn_prf <- read_sf(file.path(dir,"data","zn_prf.shp"))
#' geo_meuble <- read_sf(file.path(dir,"data","geo_meuble.shp"))
#' p.pendage <- read_sf(file.path(dir,"data","p.pendage.shp"))
#' geo_dalle(mnt0=mnt, zn_prg=zn_prf, geo_meuble=geo_meuble, p.pendage=p.pendage)
#'
geo_dalle <- function(mnt0, zn_prf=NULL, zn = NULL, geo_meuble = NULL, p.pendage = NULL,
                      con_geo = NULL, db_geo = NULL, db_pendage = NULL,
                      geo_not_dalle = c("j", "n", "c", "l"), buff = 500){


# emprise -----------------------------------------------------------------

  if(is.null(zn)){
    zn <- zn_prf %>% st_bbox() %>% st_as_sfc() %>% st_buffer(buff) %>% st_as_sf()
  }

# postgis géologie --------------------------------------------------------

if(! is.null(con_geo)){

  bb <- st_bbox(zn %>% st_transform(2154))

  geo <- st_read(con_geo, query = paste0(
    "select * from \"", db_geo,  "\" WHERE ST_Intersects (geometry, ST_MakeEnvelope (",
    bb["xmin"], ",", bb["ymin"], ",", bb["xmax"], ",", bb["ymax"], ",2154))")
  )

  non_meuble <- map(geo_not_dalle, ~which(startsWith(geo$notation,.x))) %>% unlist()
  geo_meuble <- geo %>% slice(-non_meuble)

  p.pendage <- st_read(con_geo, query = paste0(
    "select * from \"", db_pendage,  "\" WHERE ST_Intersects (geometry, ST_MakeEnvelope (",
    bb["xmin"], ",", bb["ymin"], ",", bb["xmax"], ",", bb["ymax"], ",2154))")
  )

}

  # rasters -----------------------------------------------------------------

  mnt <- mnt0 %>% crop(zn)


  if(res(mnt)[1]<5){
    mnt <- mnt %>% aggregate(fact=5/res(mnt)[1])
  }


  pte <- terrain(mnt,"slope")
  st <- stack(c(mnt= resample(mnt,pte), pte = pte))

  # rasters lissés à 20m
  fw <- focalWeight(st$mnt,20,type = "circle")
  st$mnts <- focal(st$mnt, w=fw)
  st$ptes <- terrain(st$mnts,"slope")
  st$asps <- terrain(st$mnts,"aspect")
  st$aspsx <- cos(st$asps)
  st$aspsy <- sin(st$asps)
  st$cur <- curvature(st$mnts,"bolstad")


  # suppression des pentes irrégulières
  # jugées selon curvature

  st$ind <- st$ptes
  st$ind[abs(st$cur)>.01] <- NA

  # suppression des comblements

  if(is.null(geo_meuble)){
    # jugées selon tri et pente faible

    st$tri <- tri(st$mnt, exact = FALSE)
    st$ind[((st$tri * st$pte)^.5) < .5] <- NA

  }else if(nrow(geo_meuble)>0){
    # jugées selon nature géologie (si dispo)

    st$ind <- mask(st$ind, geo_meuble %>% st_transform(st_crs(st$ind)), inverse = T)
  }

  # découpage en zones homogènes

  st$clump <- clump(st$ind)

  # scission des unités de + de 1000 px hétérogènes

  find_clu_a_reduire <- function(st){
    val_clump <- data.frame(clump=values(st$clump), pte=values(st$ptes),azx=values(st$aspsx), azy=values(st$aspsy)) %>%
      na.omit() %>% group_by(clump) %>%
      summarise(nb = n(),
                m_pte = mean(pte),
                cv_pte = sd(pte)/m_pte,
                sd_aspx=sd(azx) * (m_pte>.2),
                sd_aspy=sd(azy) * (m_pte>.2))

    val_clump %>%
      filter(nb>1000 & (cv_pte>0.1 | sd_aspx>0.05 | sd_aspy>0.05)) %>%
      pull(clump)
  }

  clu_a_reduire <- find_clu_a_reduire(st)

  while (length(clu_a_reduire)>0) {

    big_cl <- st$clump
    big_cl[! values(big_cl) %in% clu_a_reduire] <- NA
    st$clump[big_cl %>% boundaries() == 1] <- NA
    st$clump <- clump(st$clump)

    clu_a_reduire <- find_clu_a_reduire(st)
  }

  v <- values(st) %>% as.data.frame()
  v$x <- xFromCell(st$mnt, 1:ncell(st$mnt))
  v$y <- yFromCell(st$mnt, 1:ncell(st$mnt))

  # suppression des zones de moins de 50 px

  clump20 <- v %>% na.omit() %>% count(clump) %>% filter(n>50) %>% pull(clump)
  v$clump[! v$clump %in% clump20] <- NA

  # suppression des clumps avec moins de 3 valeurs de x ou y

  clump_rm <- v %>% group_by(clump) %>%
    summarise(count_x = length(unique(x)), count_y = length(unique(y))) %>%
    filter(count_x < 3 | count_y < 3) %>% pull(clump)

  v$clump[v$clump %in% clump_rm] <- NA

  # estimation du plan dans lequel se situe chaque zone

  v$clump <- v$clump %>% as.factor() %>% as.numeric()

  clu_names <- unique(v$clump) %>% na.omit()
  names(clu_names) <- clu_names


  vplan <- map(clu_names,function(i){

    m <- glm(mnt~poly(x,2)+ poly(y,2), data = v %>% na.omit() %>%
               filter(clump==i))
    p <- predict(m,v)
    p-mean(p)
  })

  vplan <- do.call(rbind, vplan) %>% as.data.frame()

  # classification des plans

  d <- dist(vplan %>% select(sample.int(ncol(vplan),1000)))
  hd <- hclust(d, method = "ward.D")
  brk_k <- classBreaks(hd$height,2,"std")[2]

  k <- sum(hd$height>=brk_k)
  cl_hd <- cutree(hd,k)

  plot(hd)
  rect.hclust(hd,k=k,border="blue")

  v$m_cl <- NULL
  v <- v %>% left_join(data.frame(clump = as.numeric(as.character(names(cl_hd))), m_cl = cl_hd))

  plan_cl <- st$mnt
  values(plan_cl) <- v$m_cl
  plotClu <- qm(plan_cl) + tm_shape(p.pendage)+tm_text("AZIMUT")

  sum_clu <- v %>% na.omit() %>%  group_by(m_cl) %>% summarise(
    nb = round(n()),
    pte_m = mean(pte),pte_sd=sd(pte),
    cur_m = mean(cur)*1000, cur_sd=sd(cur)*1000,
    mnt_m= mean(mnt), mnt_sd=sd(mnt),
    aspx = mean(aspsx), aspy = mean(aspsy)

  ) %>% arrange(desc(nb))

  # comparaison avec données pendage BRGM -----------------------------------------

  # points situés sur le même versant que la parcelle (si pente parcelle >30%)


  if(! is.null(p.pendage)){

    (val_pt_topo <- extract(st[[c("ptes","aspsx","aspsy")]],
                                    p.pendage %>% st_buffer(50),
                                    fun = mean, na.rm=TRUE) %>%
       as.data.frame() %>%
       mutate(asp = atan(aspsy/aspsx)/pi*180))
    p.pendage$asp <- val_pt_topo$asp


    val_topo_prf <- extract(st,zn_prf %>% st_union() %>% as_Spatial(), fun=mean)[1,]

    (asp_prf <- atan(val_topo_prf["aspsy"]/val_topo_prf["aspsx"])/pi*180)

    if(atan(val_topo_prf["ptes"])>0.3){
      p.pendage.v <- p.pendage %>% filter(abs(asp-asp_prf)<45) %>%
        mutate(PENDAGE = replace(PENDAGE,PENDAGE == 999,NA))
    }else{
      p.pendage.v <- p.pendage
    }

    # clump le plus prche des valeurs moyennes des pointsrestants

    if(nrow(p.pendage.v)>0){
      (diff_azx <- sum_clu$aspx - mean(cos(p.pendage.v$AZIMUT/180*pi)))
      (diff_azy <- sum_clu$aspy - mean(sin(p.pendage.v$AZIMUT/180*pi)))
      sum_clu$diff_az <- abs(diff_azx * diff_azy)

      moy_pte_pts <- mean(p.pendage.v$PENDAGE, na.rm=TRUE)
      if(!is.na(moy_pte_pts)){
        (sum_clu$diff_pte <- abs(moy_pte_pts - sum_clu$pte_m *180 / pi))
      }
    }

    # élimine les clumps très différents

    sum_clu_select <- sum_clu

    if("diff_az" %in% colnames(sum_clu)){
      sum_clu_select <- sum_clu_select %>%
        filter((diff_az/mean(diff_az))<2)
    }

    if("diff_pte" %in% colnames(sum_clu)){
      sum_clu_select <- sum_clu_select %>%
        filter((diff_pte/mean(diff_pte))<2)
    }

    choix <- sum_clu_select$m_cl[1]
  }else{
    choix <- sum_clu$m_cl[1]
  }

  # création du plan supposé de la strate géologique

  vplan_choix <- vplan[which(cl_hd %in% choix),]
  clump_choix <- names(cl_hd[cl_hd %in% choix]) %>% as.numeric()

  cnt <- v %>% na.omit() %>% count(clump) %>%
    filter(clump %in% clump_choix)

  st$plan <- st$mnt
  values(st$plan) <- apply(vplan_choix,2, weighted.mean, w = cnt$n)

  plotPlan <- qm(st$plan)

  # pente du terrain par rapport à ce plan

  st$h <- st$mnt - st$plan
  st$ph <- terrain(st$h,"slope")

  rr <-st$ph*100
  rr[rr>20] <- 20
  rr[abs(st$cur)>.02] <- 20

  # rr <- mask(rr, geo_meuble %>% st_transform(st_crs(st$ind)), inverse = T)
  rrc <- rr %>% crop(zn_prf %>% st_buffer(100) %>% as_Spatial())
  stc <- st %>% crop(zn_prf %>% st_buffer(100) %>% as_Spatial())

  plot2D <- tmap_arrange(qm(rrc) + tm_shape(zn_prf) + tm_borders(col="black"),
                         qm(stc$pte) + tm_shape(zn_prf) + tm_borders(col="black"), ncol=2)

  if(all(c("rgl","rasterVis") %in% installed.packages())){
    options(rgl.printRglwidget = TRUE) # si serveur
    rgl::rgl.clear()
    plot3D <- rasterVis::plot3D(stc$mnt,drape=rrc,col=hcl.colors, zfac=2)
  }else{
    plot3D <- NULL
    message("Pour visualiser la carte 3D, installez les packages 'rgl' et 'rasterViz'")
  }


  sum_clu$choix <- (sum_clu$m_cl == choix)

  return(list(raster = rr, plot2D = plot2D, plot3D = plot3D,
              plotPlan = plotPlan, plotClu = plotClu, plan = st$plan,
              mnt = st$mnt, pte = st$pte, sum_clu = sum_clu))
}



