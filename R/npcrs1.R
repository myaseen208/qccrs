#' @name    npcrs1
#' @aliases npcrs1
#' @title   NP Control Charts under Repetitive Sampling with single positive integer.
#' @description Calculates Average Sample Numbers (ASN), Average Run Length (ARL1) and value of k for NP control charts under repetitive sampling as given in Aslam et al.(2014)
#'
#' @param .n      Sample Size
#' @param .p0     probability that process is in control
#' @param .f      Size of the Shift
#' @param .ssize  Number of samples with replacement at each iteration
#' @param .kr     Random Positive Constant
#' @param .k      Positive Constant
#'
#'
#'
#' @return ARL0, ARL1 and K
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Muhammad Aslam  (\email{aslam_ravian@@hotmail.com})
#'          \item Sami Ullah      (\email{samiullahuos@@gmail.com})
#'          \item Muhammad Azam   (\email{mazam@@uvas.edu.pk})
#'          \item Chi-Hyuck Jun   (\email{chjun@@postech.ac.kr})
#'          \item Muhammad Kashif (\email{mkashif@@uaf.edu.pk})
#'          }
#'
#' @references
#' Aslam, M., Azam, M. and Jun, C. (2014).
#'  New Attributes and Variables Control Charts under Repetitive Sampling.
#'  \emph{Industrial Engineering & Management Systems}.
#'  \strong{13}(1):101-106.
#'
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
#' @importFrom stats dbinom runif
#' @importFrom tibble tibble
#'
#' @export
#'
#' @examples
#' 
#' library(magrittr)
#' npcrs1(
#'   .n     = 60
#' , .p0    = 0.10
#' , .f     = 0.10
#' , .k     = 2.6432
#' )
#' 
#' 
#' npcrs1(
#'   .n     = 60
#' , .p0    = 0.10
#' , .f     = 0.10
#' , .ssize = 1000
#' , .kr    = 4
#' )
#'
#' 
if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
         "ARL1"
        , "ASN"
        , "ASS0"
	, "LCL1"
	, "LCL2"
	, "LCL22"
	, "Pin"
	, "Pinc"
	, "Prep"
	, "Prep1"
	, "UCL1"
	, "UCL2"
	, "k"
	, "k1"
	, "k2"
    )
  )
}
     
npcrs1 <- function(.n, .p0, .f , .ssize = NULL, .k = NULL,  .kr = NULL){
  UseMethod("npcrs1")
}
#' @export
#' @rdname npcrs1

npcrs1.default <- function(.n, .p0, .f , .ssize = NULL, .k = NULL,  .kr = NULL){
  if(!is.null(.k) & is.null(.ssize) & is.null(.kr))
  {
    mu    <- .n*.p0
    sigma <- sqrt(.n*.p0*(1-.p0))
    .p1   <- .p0 + .f*.p0
    LCL   <- as.integer(mu - .k*sigma)
    UCL   <- as.integer(mu + .k*sigma)
    Pc    <- sum(x = dbinom((LCL+1):UCL, size = .n, prob = .p0, log = FALSE))
    Pcc   <- sum(x = dbinom((LCL+1):UCL, size = .n, prob = .p1, log = FALSE))
    ARL0  <- 1/(1-Pc)
    ARL1  <- 1/(1-Pcc)
    return(tibble::tibble(ARL0 = ARL0 , ARL1=ARL1))
    }  
  else 
    {
      mu    <- .n*.p0
      sigma <- sqrt(.n*.p0*(1-.p0))
      df1 <- 
        tibble::tibble(
          k = runif(n =  .ssize, min = 0, max = .kr)  
        ) %>%
        dplyr::mutate(
            LCL = as.integer(mu - k*sigma)
          , UCL = as.integer(mu + k*sigma)
          , LCL1 = ifelse(test = LCL < 0, yes = -1, no = LCL)
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Pc   = sum(dbinom(x = (LCL1  + 1):UCL, size = .n, prob = .p0, log = FALSE))
          , ARL0 = 1/(1-Pc)
        ) %>%
        dplyr::filter(ARL0 >= 90) %>%
        dplyr::arrange(ARL0) %>%
        dplyr::select(k, ARL0) %>% 
        dplyr::filter(dplyr::row_number()==1)
      
      .k1   <- df1$k
      ARL0  <- df1$ARL0
      mu    <- .n*.p0
      sigma <- sqrt(.n*.p0*(1-.p0))
      .p1   <- .p0 + .f*.p0
      LCL1   <- as.integer(mu - .k1*sigma)
      UCL1   <- as.integer(mu + .k1*sigma)
      Pc1    <- sum(x = dbinom((LCL1+1):UCL1, size = .n, prob = .p0, log = FALSE))
      Pcc1   <- sum(x = dbinom((LCL1+1):UCL1, size = .n, prob = .p1, log = FALSE))
      ARL00  <- 1/(1-Pc1)
      ARL10  <- 1/(1-Pcc1)
      return(tibble::tibble(k = .k1, ARL0 = ARL00 , ARL1=ARL10))
    }
  }
