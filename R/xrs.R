#' @name    xrs
#' @aliases xrs
#' @title   Xbar Control Charts Under Repetitive Sampling
#' @description Calculates the Average Sample Number and Average Run Length as given in Aslam et al. (2014)
#'
#' @param .n     Sample Size
#' @param .c     Size of the Shift      
#' @param .k1    Positive Integer
#' @param .k2    Positive Integer
#'
#'
#' @return  Average Sample Number (ASN) and Average Run Length (ARL1) for xbar control charts under repetitive sampling
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Muhammad Aslam (\email{aslam_ravian@@hotmail.com})
#'          \item Sami Ullah (\email{samiullahuos@@gmail.com})
#'          \item Muhammad Azam   (\email{mazam@@uvas.edu.pk})
#'          \item Chi-Hyuck Jun   (\email{chjun@@postech.ac.kr})
#'          \item Muhammad Kashif (\email{mkashif@@uaf.edu.pk})
#'          }
#'
#' @references
#' Aslam, M., Azam, M. and Jun, C. (2014).
#'  New Attributes and Variables Control Charts under Repetitive Sampling.
#'  \emph{Industrial Engineering & Management Systems}.
#'  \strong{1}(13):101-106.
#'
#' @importFrom stats pnorm 
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map
#'
#' @export
#'
#' @examples
#' 
#' library(magrittr)
#' library(purrr)
#' 
#' c(0.0, 0.1, 0.20, 0.3, 0.4, 0.5, 1.0, 1.5, 2, 3) %>%
#' purrr::map(
#' function(x) 
#'     xrs(
#'         .c     = x
#'       , .n     = 10
#'       , .k1    = 2.9301
#'       , .k2    = 0.9825))
#'
xrs <- function(.c, .n,  .k1 , .k2 ){
  UseMethod("xrs")
}

#' @export
#' @rdname xrs
xrs.default <- function(.c, .n,  .k1 , .k2 ){
  df <-  tibble::tibble(
       .c     = .c  
      , Pin   = (2* pnorm(q = .k2, mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-1)/
              1- 2*(pnorm(q = .k1, mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-
                    pnorm(q = .k2, mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE))
      , Pinc  = (pnorm(q = .k2 - .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
                 pnorm(q = .k2 + .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-1)/
                (pnorm(q = .k2 + .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-
                 pnorm(q = .k1 + .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-
                 pnorm(q = .k1 - .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
                 pnorm(q = .k2 - .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)+1)
      , Prep  = 2*(pnorm(q = .k1, mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-
                   pnorm(q = .k2, mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE))
      , Prep1 = -pnorm(q = .k2 + .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
                 pnorm(q = .k1 + .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
                 pnorm(q = .k1 - .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)-
                 pnorm(q = .k2 - .c*sqrt(.n), mean =0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      , ARL0  =   1/(1-Pin)
      , ARL1  =   1/(1-Pinc)
      , ASN  =   .n/(1-Prep1)
    ) %>%
      dplyr::select(.c, ASN, ARL1)
    return(df)
  }
