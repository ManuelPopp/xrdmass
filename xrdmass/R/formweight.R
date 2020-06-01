#' form.weight
#'
#' Calculates the molar weight of a chemical compound or the sum of the atomic weights for each element of a compound.
#' @param formula A character string containing a chemical formula, e.g. "SiO2" or "Ca(OH)2"
#' @param out Defines the output format of the function. If unspecified or out = "formula" the function gives out the molecular weight for the given formula. If out = "element" the output is a data.frame which lists all elements and the respective weight multiplied by the number of atoms of this element within the compound.
#' @param update Binary (TRUE or FALSE). Sets whether the information on the elements and their atomic masses shall be updated from GitHub. The option might be outdated and the information not longer available online; if not set, the function uses an internal database.
#' @return A numeric value or data.frame, as defined by the function parameter "out".
#' @export
form.weight <- function(formula, out = "formula", update = FALSE){

  ################################################
  #### handle brackets
  ################################################
  # check for typos
  if(NROW(unlist(gregexpr("[\\(]", formula))) != NROW(unlist(gregexpr("[\\)]", formula)))){
    if(NROW(unlist(gregexpr("[\\(]", formula))) > NROW(unlist(gregexpr("[\\)]", formula)))){
      stop("Formula contains unmatched opening bracket.")
    }else if(NROW(unlist(gregexpr("[\\(]", formula))) < NROW(unlist(gregexpr("[\\)]", formula)))){
      stop("Formula contains unmatched closing bracket.")
    }
  }
  caps <- unlist(gregexpr("[[:upper:]]", formula))

  while(unlist(gregexpr("[[:upper:]][[:upper:]]", formula))[1] >= 1){
    up <- unlist(gregexpr("[[:upper:]][[:upper:]]", formula))
    frml <- paste(
      substr(formula, 1, up),
      "1",
      substr(formula, up+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  while(unlist(gregexpr("[[:lower:]][[:upper:]]", formula))[1] >= 1){
    low <- unlist(gregexpr("[[:lower:]][[:upper:]]", formula))
    frml <- paste(
      substr(formula, 1, low),
      "1",
      substr(formula, low+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  while(unlist(gregexpr("[[:upper:]]\\(", formula))[1] >= 1){
    brac1 <- unlist(gregexpr("[[:upper:]]\\(", formula))
    frml <- paste(
      substr(formula, 1, brac1),
      "1",
      substr(formula, brac1+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  while(unlist(gregexpr("[[:lower:]]\\(", formula))[1] >= 1){
    brac2 <- unlist(gregexpr("[[:lower:]]\\(", formula))
    frml <- paste(
      substr(formula, 1, brac2),
      "1",
      substr(formula, brac2+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  while(unlist(gregexpr("[[:upper:]]\\)", formula))[1] >= 1){
    brac3 <- unlist(gregexpr("[[:upper:]]\\)", formula))
    frml <- paste(
      substr(formula, 1, brac3),
      "1",
      substr(formula, brac3+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  while(unlist(gregexpr("[[:lower:]]\\)", formula))[1] >= 1){
    brac4 <- unlist(gregexpr("[[:lower:]]\\)", formula))
    frml <- paste(
      substr(formula, 1, brac4),
      "1",
      substr(formula, brac4+1, nchar(formula)),
      sep = ""
    )
    formula <- frml
  }

  if(!grepl("\\d", substr(formula, nchar(formula), nchar(formula)))){
    formula <- paste(formula, "1", sep = "")
  }

  if(grepl("(", formula, fixed = TRUE)){
    # list pairs of matching brackets
    pairs <- list()
    form <- formula
    opening_brackets <- unlist(gregexpr("[\\(]", form))
    for(i in 1:NROW(opening_brackets)){
      opening_brackets <- unlist(gregexpr("[\\(]", form))
      closing_brackets <- unlist(gregexpr("[\\)]", form))
      brackets <- unlist(
        regmatches(form, gregexpr("\\b(\\(|\\))\\b", form))
      )
      open <- unlist(
        gregexpr("\\(\\)", paste(brackets, collapse = ""))
      )[1]
      pair <- c(
        opening_brackets[open],
        closing_brackets[
          which(closing_brackets > opening_brackets[open])[1]
          ]
      )
      pairs[[i]] <- pair
      substr(form, pair[1], pair[1]) <- " "
      substr(form, pair[2], pair[2]) <- " "
    }

    # get multiplication factors for brackets
    factors <- list()
    for(i in 1:NROW(pairs)){
      start <- pairs[[i]][2]+1
      x <- unlist(gregexpr("[\\(|\\)|[:upper:]]", formula))
      if(start < rev(x)[1]){
        end <- x[which(x > start)[1]]-1
      }else{
        end <- nchar(formula)
      }
      if(!is.na(as.numeric(substr(formula, start, end)))){
        factors[[i]] <- substr(formula, start, end)
        substr(formula, start, end+1) <- paste(rep(" ", end+1-start), collapse = "")
      }else{
        pairs <- pairs[-i]
      }
    }

    # remove unwanted values
    if(length(which(is.na(factors))) > 0){
      factors[[which(is.na(factors))]] <- "1"
    }
    if(length(which(factors == "")) > 0){
      factors[[which(factors == "")]] <- "1"
    }

    options(digits = 9)
    factors <- as.numeric(factors)

    # multiply brackets
    strt <- unlist(gregexpr("[[:digit:]|\\.]+", formula))
    nd <- strt + unlist(attr(gregexpr("[[:digit:]|\\.]+", formula)[[1]],"match.length"))-1
    vls <- unlist(regmatches(formula, gregexpr("[[:digit:]|\\.]+", formula)))
    df <- data.frame(start = strt, end = nd, value = vls)
    options(digits = 9)
    df[, 3] <- as.character(df[, 3])
    df[, 3] <- as.numeric(df[, 3])

    for(i in 1:NROW(pairs)){
      df[which(
        df[, 1] > pairs[[i]][1] & df[, 2] < pairs[[i]][2]
      ), 3] <- df[which(
        df[, 1] > pairs[[i]][1] & df[, 2] < pairs[[i]][2]
      ), 3]*factors[[i]]
    }

    for(i in 1:NROW(vls)){
      formX <- paste(
        substr(formula, 1, unlist(gregexpr("[[:digit:]|\\.]+", formula))[i]-1),
        as.character(df[i, 3]),
        substr(formula,
               unlist(gregexpr("[[:digit:]|\\.]+", formula))[i] + unlist(
                 attr(gregexpr("[[:digit:]|\\.]+", formula)[[1]],"match.length")
               )[i], nchar(formula)),
        sep = "")
      formula <- formX
    }

    formula <- gsub(" ", "", formula)
    formula <- gsub("\\(", "", formula)
    formula <- gsub("\\)", "", formula)
  }

  ################################################
  #### generate matrix
  ################################################
  # seperate elements
  breaks <- unlist(gregexpr("[[:upper:]]", formula))
  list <- list()
  for(i in 1:NROW(breaks)){
    if(i < NROW(breaks)){
      list[[i]] <- substring(formula, breaks[i], breaks[i+1]-1)
    }else{
      list[[i]] <- substring(formula, breaks[i], nchar(formula))
    }
  }

  dfx <- data.frame(data = unlist(list))
  dfx[, 1] <- as.character(dfx[, 1])
  dfx$Element <- sub("^([[:alpha:]]*).*", "\\1", dfx[, 1])
  dfx$Count <- as.numeric(unlist(regmatches(dfx[, 1],
                                            gregexpr("[[:digit:]]+\\.*[[:digit:]]*", dfx[, 1]))))

  ################################################
  #### get element information
  ################################################
  if(update == TRUE){
    #options(warn = -1)
    #if(!require("RCurl")){
    #  install.packages("RCurl")
    #}
    #requireNamespace("RCurl")
    #tmp <- getURL("https://gist.githubusercontent.com/binary-bisam/3f92c2f77729541ec7b626197c7c6dc4/raw/1d92663004489a5b6926e944c1b3d9ec5c40900e/Periodic%2520Table%2520of%2520Elements.csv")
    #el_info <- read.table(text = tmp, sep = ",", header = TRUE)
    #element_info <- el_info[,c(1, 2, 3, 4)]
  }else{
    element_info <- data.frame(AtomicNumber = 1:74)
    element_info$Element <- structure(c(28L, 26L, 34L, 6L, 7L, 11L, 43L, 44L, 21L, 40L, 58L,
                                        36L, 1L, 56L, 46L, 60L, 14L, 3L, 47L, 10L, 54L, 67L, 68L, 15L,
                                        37L, 31L, 16L, 41L, 17L, 73L, 23L, 24L, 4L, 55L, 8L, 32L, 51L,
                                        59L, 72L, 74L, 42L, 38L, 62L, 52L, 50L, 45L, 57L, 9L, 29L, 66L,
                                        2L, 63L, 30L, 70L, 13L, 5L, 33L, 12L, 48L, 39L, 49L, 53L, 20L,
                                        22L, 64L, 18L, 27L, 19L, 65L, 71L, 35L, 25L, 61L, 69L), .Label = c("Aluminum",
                                                                                                           "Antimony", "Argon", "Arsenic", "Barium", "Beryllium", "Boron",
                                                                                                           "Bromine", "Cadmium", "Calcium", "Carbon", "Cerium", "Cesium",
                                                                                                           "Chlorine", "Chromium", "Cobalt", "Copper", "Dysprosium", "Erbium",
                                                                                                           "Europium", "Fluorine", "Gadolinium", "Gallium", "Germanium",
                                                                                                           "Hafnium", "Helium", "Holmium", "Hydrogen", "Indium", "Iodine",
                                                                                                           "Iron", "Krypton", "Lanthanum", "Lithium", "Lutetium", "Magnesium",
                                                                                                           "Manganese", "Molybdenum", "Neodymium", "Neon", "Nickel", "Niobium",
                                                                                                           "Nitrogen", "Oxygen", "Palladium", "Phosphorus", "Potassium",
                                                                                                           "Praseodymium", "Promethium", "Rhodium", "Rubidium", "Ruthenium",
                                                                                                           "Samarium", "Scandium", "Selenium", "Silicon", "Silver", "Sodium",
                                                                                                           "Strontium", "Sulfur", "Tantalum", "Technetium", "Tellurium",
                                                                                                           "Terbium", "Thulium", "Tin", "Titanium", "Vanadium", "Wolfram",
                                                                                                           "Xenon", "Ytterbium", "Yttrium", "Zinc", "Zirconium"), class = "factor")
    element_info$Symbol <- structure(c(26L, 27L, 35L, 7L, 5L, 9L, 40L, 46L, 21L, 44L, 41L,
                                       37L, 2L, 58L, 47L, 54L, 13L, 3L, 32L, 10L, 56L, 66L, 68L, 15L,
                                       38L, 22L, 14L, 45L, 17L, 73L, 23L, 25L, 4L, 57L, 8L, 33L, 51L,
                                       61L, 71L, 74L, 42L, 39L, 64L, 53L, 52L, 48L, 1L, 11L, 31L, 60L,
                                       55L, 65L, 30L, 70L, 16L, 6L, 34L, 12L, 50L, 43L, 49L, 59L, 20L,
                                       24L, 63L, 18L, 29L, 19L, 67L, 72L, 36L, 28L, 62L, 69L), .Label = c("Ag",
                                                                                                          "Al", "Ar", "As", "B", "Ba", "Be", "Br", "C", "Ca", "Cd", "Ce",
                                                                                                          "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", "Fe", "Ga",
                                                                                                          "Gd", "Ge", "H", "He", "Hf", "Ho", "I", "In", "K", "Kr", "La",
                                                                                                          "Li", "Lu", "Mg", "Mn", "Mo", "N", "Na", "Nb", "Nd", "Ne", "Ni",
                                                                                                          "O", "P", "Pd", "Pm", "Pr", "Rb", "Rh", "Ru", "S", "Sb", "Sc",
                                                                                                          "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Ti", "Tm",
                                                                                                          "V", "W", "Xe", "Y", "Yb", "Zn", "Zr"), class = "factor")
    element_info$AtomicMass <- c(1.007, 4.002, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999,
                                 18.998, 20.18, 22.99, 24.305, 26.982, 28.086, 30.974, 32.065,
                                 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996,
                                 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.64,
                                 74.922, 78.96, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
                                 92.906, 95.96, 98, 101.07, 102.906, 106.42, 107.868, 112.411,
                                 114.818, 118.71, 121.76, 127.6, 126.904, 131.293, 132.905, 137.327,
                                 138.905, 140.116, 140.908, 144.242, 145, 150.36, 151.964, 157.25,
                                 158.925, 162.5, 164.93, 167.259, 168.934, 173.054, 174.967, 178.49,
                                 180.948, 183.84)
  }

  ################################################
  #### add weights and sum up
  ################################################
  for(i in 1:nrow(dfx)){
    if(!dfx$Element[i] %in% element_info$Symbol){
      stop(paste(
        "Invalid element symbol '", dfx$Element[i], "'.", sep = ""
      ))
    }
  }

  dfx$Weight <- rep(-999, nrow(dfx))
  for(i in 1:nrow(dfx)){
    dfx$Weight[i] <- element_info[which(element_info$Symbol == dfx$Element[i]),4]
  }

  dfx$rowsum <- as.numeric(dfx$Count)*as.numeric(dfx$Weight)
  sum <- sum(dfx$rowsum)
  if(out == "formula"|out == "f"){
    return(sum)
  }
  if(out == "element"|out == "e"){
    elements <- unique(dfx$Element)
    tbl <- data.frame(element = elements)
    tbl$mass <- unique(dfx$Weight)
    tbl$count <- rep(NA, nrow(tbl))
    for(i in 1:nrow(tbl)){
      tbl$count[i] <- sum(dfx$Count[which(dfx$Element == tbl$element[i])])
    }
    tbl$sum <- tbl$mass*tbl$count
  }
  return(tbl)
}
