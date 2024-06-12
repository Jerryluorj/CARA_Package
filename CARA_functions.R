library(doParallel)
library(foreach)

###############################################################################
####################   Zhao's NEW CARA Design   ###############################
###############################################################################
#' Title
#'
#' @param PredictiveCov a vector. It specifies patients' predictive covariate \code{X} and must be discretized
#' @param PrognosticCov a matrix. It specifies patients' \code{k} prognostic covariates code{Z} and must be discretized
#' @param mu a vector of length 2. It specifies the predetermined model coefficients of treatment effect \code{T} and \code{1-T} to generate response.
#' @param theta a vector of length 2. It specifies the predetermined model coefficients of \code{X} and \code{X*T} to generate response.
#' @param gamma a vector of length \code{k}. It specifies the predetermined model coefficients of \code{Z} and \code{X*T} to generate response.
#' @param m0 an integer. If specifies the initial number of \code{2m0} patients that are assigned to treatments 1 and 2 with restrict randomization.
#' @param p an integer between 0.75 and 0.95. It is the probability of Efron's biased coin to assign patients to treatments 1 and 2.
#' @param omega a vector of length \code{k+2} that sums to 1. It specifies the weights placed on overall, on \code{k} covariate margin and within a stratum cell, respectively.
#' @param isResponseBinary a logical evaluating to TRUE or FALSE indicating whether the response is binary.
#' @param isTreatmentNumber a logical evaluating to TRUE or FALSE indicating whether the Treatment is denoted as 1 and 0.
#'
#' @return \item{data}{The patients' data after assignment.}
#' @export
#'
#' @examples ZhaoNew(PredictiveCov=sample(c(1, -1), 1000, replace = T),PrognosticCov=replicate(3, sample(c(1, -1), 1000, replace = T)),
#'  beta=c(0.5,0.5,1,-0.5, 0.5, 0.5,0.5), isResponseBinary = T, isTreatmentNumber=T, omega=rep(0.2,5))

ZhaoNew = function(PredictiveCov ,
                   PrognosticCov,
                   mu,
                   theta ,
                   gamma ,
                   m0 = 10,
                   p = 0.8,
                   omega,
                   isResponseBinary ,
                   isTreatmentNumber) {
  n = length(PredictiveCov)
  k = ncol(PrognosticCov)
  t = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  beta=c(mu,theta,gamma)
  if (isResponseBinary) {
    py = 1 / (1 + exp(
      -(
        beta[1] * t + beta[2] * (1 - t) + beta[3] * PredictiveCov + beta[4] * PredictiveCov * t +
          rowSums(beta[5:length(beta)] * PrognosticCov)
      )
    ))
    py_2m0 = py[1:2 * m0]
    Y = c(ifelse(runif(2 * m0) < py_2m0, 1, 0), rep(NA, n - 2 * m0))
  } else{
    Y = c((
      beta[1] * t + beta[2] * (1 - t) + beta[3] * PredictiveCov + beta[4] * PredictiveCov * t +
        rowSums(beta[5:length(beta)] * PrognosticCov) + rnorm(n)
    )[1:(2 * m0)],
    rep(NA, n - 2 * m0)
    )
  }
  data = cbind(PredictiveCov, PrognosticCov, Y, t)
  colnames(data)[1] = "X"
  colnames(data)[2:(1 + k)] = paste0("Z", 1:k)
  colnames(data)[2 + k] = "Y"
  colnames(data)[3 + k] = "t"
  for (r in (2 * m0 + 1):nrow(data)) {
    Xn = data[r, 1]
    Zn = data[r, 2:(1 + k)]
    data_before = data[1:(r - 1), ]

    if (isResponseBinary) {
      p1hat = sum(data_before[, "Y"] == 1 &
                    data_before[, "X"] == Xn &
                    data_before[, "t"] == 1) / sum(data_before[, "X"] == Xn &
                                                     data_before[, "t"] == 1)
      p0hat = sum(data_before[, "Y"] == 1 &
                    data_before[, "X"] == Xn &
                    data_before[, "t"] == 0) / sum(data_before[, "X"] == Xn &
                                                     data_before[, "t"] == 0)

      rho = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p0hat))
    } else{
      data_before_Xn = data_before[data_before[, "X"] == Xn, ]

      lm_formula = formula(paste0("Y ~ t + X + t:X+", paste(colnames(data)[2:(1 +
                                                                                k)], collapse = "+")))
      model = lm(lm_formula, data = data.frame(data_before_Xn))
      coefficients = coef(model)

      mu1_hat = coefficients['t']
      intercept_hat = coefficients['(Intercept)']
      beta1_hat = coefficients['X']
      beta2_hat = coefficients['t:X']
      gamma_hat = numeric()
      for (i in 1:k) {
        gamma_hat = append(gamma_hat, coefficients[paste0("Z", k)])
      }

      Yhatn_1 = intercept_hat + mu1_hat + beta1_hat * Xn + beta2_hat * Xn + sum(gamma_hat *
                                                                                  Zn)
      Yhatn_0 = intercept_hat + beta1_hat * Xn + sum(gamma_hat * Zn)

      rho = pnorm(Yhatn_1) / (pnorm(Yhatn_1) + pnorm(Yhatn_0))
    }

    data_before_1 = rbind(data_before, data[r, ])
    data_before_1[r, "t"] = 1
    data_before_2 = rbind(data_before, data[r, ])
    data_before_2[r, "t"] = 0

    Dn1 = sum(data_before_1[, "t"] == 1 &
                data_before_1[, "X"] == Xn) - rho * sum(data_before_1[, "X"] == Xn)
    Dn2 = sum(data_before_2[, "t"] == 1 &
                data_before_2[, "X"] == Xn) - rho * sum(data_before_2[, "X"] == Xn)

    Dn1k = numeric()
    Dn2k = numeric()
    for (cov in 1:k) {
      append(
        Dn1k,
        sum(
          data_before_1[, "t"] == 1 &
            data_before_1[, "X"] == Xn &
            data_before_1[, paste0("Z", cov)] == Zn[cov]
        ) -
          rho * sum(data_before_1[, "X"] == Xn &
                      data_before_1[, paste0("Z", cov)] == Zn[cov])
      )

      append(
        Dn2k,
        sum(
          data_before_2[, "t"] == 1 &
            data_before_2[, "X"] == Xn &
            data_before_2[, paste0("Z", cov)] == Zn[cov]
        ) -
          rho * sum(data_before_2[, "X"] == Xn &
                      data_before_2[, paste0("Z", cov)] == Zn[cov])
      )
    }

    Dn1k1k2 = sum(data_before_1[, "t"] == 1 &
                    data_before_1[, "X"] == Xn &
                    setequal(data_before_1[, paste0("Z", 1:k)], Zn)) -
      rho * sum(data_before_1[, "X"] == Xn &
                  setequal(data_before_1[, paste0("Z", 1:k)], Zn))

    Dn2k1k2 = sum(data_before_2[, "t"] == 1 &
                    data_before_2[, "X"] == Xn &
                    setequal(data_before_2[, paste0("Z", 1:k)], Zn)) -
      rho * sum(data_before_2[, "X"] == Xn &
                  setequal(data_before_2[, paste0("Z", 1:k)], Zn))

    Imbn1 = omega[1] * Dn1 ^ 2 + sum(omega[2:(length(omega) - 1)] * Dn1k ^
                                       2) + omega[length(omega)] * Dn1k1k2 ^ 2
    Imbn2 = omega[1] * Dn2 ^ 2 + sum(omega[2:(length(omega) - 1)] * Dn2k ^
                                       2) + omega[length(omega)] * Dn2k1k2 ^ 2

    if (is.na(Imbn1 > Imbn2) |
        is.na(Imbn1 > Imbn2) |
        is.na(Imbn1 == Imbn2)) {
      data[r, "t"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 > Imbn2) {
      data[r, "t"] = sample(c(1, 0), 1, prob = c(1 - p, p))
    }
    else if (Imbn1 == Imbn2) {
      data[r, "t"] = sample(c(1, 0), 1, prob = c(0.5, 0.5))
    }
    else if (Imbn1 < Imbn2) {
      data[r, "t"] = sample(c(1, 0), 1, prob = c(p, 1 - p))
    }


    # py=1/(1+exp(-(0.5+1*data[r,1]-0.5*data[r,1]*data[r,5]+0.5*data[r,2]+0.5*data[r,3])))
    if (isResponseBinary) {
      py = 1 / (1 + exp(-(
        beta[1] * data[r, "t"] + beta[2] * (1 - data[r, "t"]) + beta[3] * Xn + beta[4] * Xn * data[r, "t"] +
          sum(beta[5:length(beta)] * Zn)
      )))
      data[r, "Y"] = sample(c(1, 0), 1, prob = c(py, 1 - py))
    } else {
      data[r, "Y"] = beta[1] * data[r, "t"] + beta[2] * (1 - data[r, "t"]) + beta[3] * Xn + beta[4] * Xn * data[r, "t"] +
        sum(beta[5:length(beta)] * Zn) + rnorm(1)
    }

  }
  if (isTreatmentNumber) {
    return(data)
  } else {
    data[, "t"] = ifelse(data[, "t"] == 1, "A", "B")
    return(data.frame(data))
  }
}
