#' XRD mass
#'
#' Calculates the mass for the elements in a sample based on a data.frame which contains the compound formulas
#' and the respective percentual share of each compound.
#' @param data A data.frame which contains the chemical compound formulas in the first column and the respective share in the second column.
#' @return A data.frame which contains the elements found in the sample and their percentual share by mass.
#' @export
################################################
#### calculate percentual mass of each element
#### based on formula and percentual share of the
#### respective phase
################################################
xmass <- function(data){
  if(sum(data[,2] <= 1)){
    p_format <- 1
  }else if(sum(data[,2] <= 100)){
    p_format <- 100
  }else{
    warning("Sum of phases > 100%.")
    p_format <- 100
  }
  data[,1] <- as.character(data[,1])
  lst <- list()
  for(i in 1:nrow(data)){
    fweight <- form.weight(formula = data[i,1], out = "element")
    pweight <- fweight
    pweight$psum <- fweight$sum*data[i,2]/p_format
    lst[[i]] <- pweight
  }

  all_phases <- do.call("rbind",lst)
  all_phases$element <- as.character(all_phases$element)
  out <- data.frame(element = unique(all_phases$element))
  out$element <- as.character(out$element)
  out$mass <- unique(all_phases$mass)
  out$percentual_mass <- rep(NA, nrow(out))
  for(i in 1:nrow(out)){
    out$percentual_mass[i] <- sum(all_phases$psum[which(all_phases$element == out$element[i])])
  }
  return(out)
}
