#' @name    npcrs2
#' @aliases npcrs2
#' @title   Attributes Control Charts under Repetitive Sampling with two positive integers
#' @description Calculates Average Sample Numbers (ASN), Average Run Length (ARL1) and value of k1 and k2 for attributes control charts under repetitive sampling as given in Aslam et al.(2014)
#'
#' @param .p0     probability that process is in control
#' @param .f      Size of the Shift
#' @param .n      Sample Size
#' @param .ssize  Number of samples with replacement at each iteration
#' @param .k1r    Random postive constant
#' @param .k2r    Random postive constant
#' @param .k1     Fixed positive constant
#' @param .k2     Fixed positive constant
#'
#'
#'
#' @return ASN, ARL, K1 and K2
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
#'  npcrs2(
#'   .n     = 40
#' , .p0    = 0.10
#' , .f     = 0.1
#' , .ssize = 1000
#' , .k1r   = 4
#' , .k2r   = .95
#'    )
#'
#'
#' npcrs2(
#'   .n     = 40
#' , .p0    = 0.10
#' , .f     = 0.1
#' , .k1    = 3.13
#' , .k2    = .731
#'    )
#'
#'
#'
#'
npcrs2 <- function(.n, .p0, .f , .ssize = NULL, .k1 = NULL, .k2 = NULL, .k1r = NULL, .k2r = NULL){
  UseMethod("npcrs2")
}

#' @export
#' @rdname npcrs2



npcrs2.default <- function(.n, .p0, .f , .ssize = NULL, .k1 = NULL, .k2 = NULL, .k1r = NULL, .k2r = NULL){
  if(!is.null(.k1) & !is.null(.k2) & is.null(.ssize) & is.null(.k1r) & is.null(.k2r))
  {
    df2 <-  tibble::tibble(
      .f    = .f  
    ,.p1    = .p0 + .f*.p0
    , mu    = .n*.p0
    , sigma = sqrt(.n*.p0*(1-.p0))
    , LCL1  = as.integer(mu - .k1*sigma)
    , UCL1  = as.integer(mu + .k1*sigma)
    , LCL2  = as.integer(mu - .k2*sigma)
    , UCL2  = as.integer(mu + .k2*sigma)
    , Pc    = sum(dbinom(x = (LCL2+1):(UCL2), size = .n, prob = .p0, log = FALSE))
    , Pcc   = sum(dbinom(x = (LCL2+1):(UCL2), size = .n, prob = .p1, log = FALSE))
    , Prep0 = sum(dbinom(x = (LCL1+1):(LCL2), size = .n, prob = .p0, log = FALSE)) + 
              sum(dbinom(x = (UCL2+1):(UCL1), size = .n, prob = .p0, log = FALSE))
    , Prep1 = sum(dbinom(x = (LCL1+1):(LCL2), size = .n, prob = .p1, log = FALSE)) +
              sum(dbinom(x = (UCL2+1):(UCL1), size = .n, prob = .p1, log = FALSE))
    , Pin0  = Pc/(1-Prep0)
    , Pin1  = Pcc/(1-Prep1)
    , ARL0  =   1/(1-Pin0)
    , ARL1  =   1/(1-Pin1)
    , ASN0  =   .n/(1-Prep0)
    , ASN1  =   .n/(1-Prep1)
    ) %>%
      dplyr::select(.f, ASN1, ARL1)
    return(df2)
  }
  else 
    {
     df1 <- 
  tibble::tibble(
    mu    =  .n*.p0
  , sigma =  sqrt(.n*.p0*(1-.p0))
  , k1    = runif(n =  .ssize, min = 3, max = .k1r)
  , k2    = runif(n =  .ssize, min = 0, max = .k2r)  
  ) %>%
  dplyr::filter(k1 > k2) %>%
  dplyr::mutate(
      LCL1 = mu - k1*sigma
    , UCL1 = mu + k1*sigma
    , LCL2 = mu - k2*sigma
    , UCL2 = mu + k2*sigma
    , LCL11 = dplyr::if_else(condition = LCL1 < 0, true = -1, false = LCL1, missing = NULL)
    , LCL22 = dplyr::if_else(condition = LCL2 < 0, true = -1, false = LCL2, missing = NULL)  
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
      Pc   = sum(dbinom(x = as.integer((LCL22 + 1)):as.integer(UCL2 ), size = .n, prob = .p0, log = FALSE))
    , Prep = sum(dbinom(x = as.integer((LCL11 + 1)):as.integer(LCL22), size = .n, prob = .p0, log = FALSE)) + 
             sum(dbinom(x = as.integer((UCL2  + 1)):as.integer(UCL1 ), size = .n, prob = .p0, log = FALSE))
    , Pin0 = Pc/(1 - Prep)
    , ARL0 = 1/(1-Pin0)
    , ASS0 =  .n/(1 - Prep)
  ) %>%
  dplyr::filter(min(ASS0) & ARL0 >= 90) %>%
  dplyr::select(k1, k2, ARL0) %>%
  dplyr::arrange(ARL0) %>% 
  dplyr::filter(dplyr::row_number()==1)
  .k1   <- df1$k1
  .k2   <- df1$k2
  .ARL0 <- df1$ARL0
  mu    <-  .n*.p0
  sigma <-  sqrt(.n*.p0*(1-.p0))
  .f     <- .f
  .p1    <- .p0 + .f*.p0
  LCL11  <- as.integer(mu - .k1*sigma)
  UCL11  <- as.integer(mu + .k1*sigma)
  LCL21  <- as.integer(mu - .k2*sigma)
  UCL21  <- as.integer(mu + .k2*sigma)
  Pc     <- sum(dbinom(x = (LCL21+1):(UCL21), size = .n, prob = .p0, log = FALSE))
  Pcc    <- sum(dbinom(x = (LCL21+1):(UCL21), size = .n, prob = .p1, log = FALSE))
  Prep0  <- sum(dbinom(x = (LCL11+1):(LCL21), size = .n, prob = .p0, log = FALSE)) + 
          sum(dbinom(x = (UCL21+1):(UCL11), size = .n, prob = .p0, log = FALSE))
  Prep1 <- sum(dbinom(x = (LCL11+1):(LCL21), size = .n, prob = .p1, log = FALSE)) +
          sum(dbinom(x = (UCL21+1):(UCL11), size = .n, prob = .p1, log = FALSE))
  Pin0  <- Pc/(1-Prep0)
  Pin1  <- Pcc/(1-Prep1)
  ARL0  <-   1/(1-Pin0)
  ARL1  <-   1/(1-Pin1)
  ASN0  <-   .n/(1-Prep0)
  ASN1  <-   .n/(1-Prep1)
  df2   <- tibble::tibble(.k1, .k2, .f, ASN1, ARL1)
  return(df2)
  }
}
