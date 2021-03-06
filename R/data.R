#' Fungal Pathogen Data
#'
#' A data set containing detection of a fungal pathogen's DNA measured at each
#' of 20 ponds located in Coconino National Forest, Arizona, USA. At each site,
#' water samples were collected to extract the DNA of Batrachochytrium
#' dendrobatidis (Bd), a fungal pathogen of various species of amphibians.
#' Environmental and ecological variables were also measured at each of the 20
#' ponds, and are included in this dataset.
#'
#' @format A data frame with 80 rows and 8 variables:
#' \describe{
#'   \item{site}{location sampled}
#'   \item{sample}{sample number from site}
#'   \item{pcr1}{detection (1) or non-dection (0) from first PCR replicate}
#'   \item{pcr2}{detection (1) or non-dection (0) from second PCR replicate}
#'   \item{loadT1}{average fungal pathogen load of frogs that tested positive for Bd and were captured during 1st survey}
#'   \item{loadT2}{average fungal pathogen load of frogs that tested positive for Bd and were captured during 2nd survey}
#'   \item{frogs}{index of frog density (catch per unit effort)}
#'   \item{Bd}{prevalence of Bd among frogs captured during 1st or 2nd surveys}
#'   }
#' @source Schmidt BR, Kery M, Ursenbacher S, Hyman OJ, and Collins JP (2013)
#'   Site occupancy models in the analysis of environmental DNA presence/absence
#'   surveys: a case study of an emerging amphibian pathogen. Methods in Ecology
#'   and Evolution 4: 646-653.
#'
"fung"

#' Tidewater Goby Data
#'
#' A data set containing detection of tidewater goby DNA in water samples
#' collected at each of 39 sites located along 400 km of coastline in California
#' and Oregon, USA. At each site, water samples were collected to extract the
#' DNA of the tidewater goby (Eucyclogbobius newberryi), a federally endangered
#' fish species living in estuaries, lagoons, and sloughs. Environmental and
#' ecological variables were also measured at each of the 39 sites, and are
#' included in this dataset.
#'
#' @format A data frame with 356 rows and 13 variables: \describe{
#'   \item{site}{location sampled}
#'   \item{sample}{sample number from site}
#'   \item{pcr1}{detection (1) or non-dection (0) from first PCR replicate}
#'   \item{pcr2}{detection (1) or non-dection (0) from second PCR replicate}
#'   \item{pcr3}{detection (1) or non-dection (0) from third PCR replicate}
#'   \item{pcr4}{detection (1) or non-dection (0) from fourth PCR replicate}
#'   \item{pcr5}{detection (1) or non-dection (0) from fifth PCR replicate}
#'   \item{pcr6}{detection (1) or non-dection (0) from sixth PCR replicate}
#'   \item{twg}{abundance index of tidewater gobies (catch per unit effort)}
#'   \item{sal}{salinity (parts per thousand)}
#'   \item{turb}{turbidity (measured as water filtration time)}
#'   \item{fish}{abundance index of fishes other than tidewater gobies (catch
#'   per unit effort)}
#'   \item{veg}{binary indicator of presence (Pres) or absence (Abs) of
#'   vegetation (widgeongrass and filamentous algae)}
#'   }
#' @source Schmelzle MC, Kinziger AP (2015) Data from: Using occupancy modeling
#'   to compare environmental DNA to traditional field methods for
#'   regional-scale monitoring of an endangered aquatic species.
"goby"
