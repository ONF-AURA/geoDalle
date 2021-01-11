

# cartes rapides

#' Carte d'un raster
#'
#' @param r raster
#' @param n nb de couleurs
#' @param pal palette de couleurs
#' @param title titre
#' @param breaks vecteur des limites de couleurs
#'
#' @return carte tmap
#'
#' @importFrom tmap tm_shape tm_raster tm_layout
#' @export
#'
qm <- function(r, n=NULL, pal="Spectral", title=NULL, breaks=NULL){
  nval <- length(unique(values(r)))

  if(!is.null(breaks)){

  }else{
    if(is.null(n)){
      if(nval<10) n <- nval else n <- 10
    }
  }
  tm_shape(r)+tm_raster(palette=pal,n=n,breaks = breaks)+
    # tm_shape(zn_prf)+tm_borders(col="black") +
    tm_layout(main.title = title)
}

