

#permuted block assignment for first 2m0 patients
permuted_block = function(n, block_size = 4) {
    assignments = numeric(n)
    treatments = c(1, 0)

    num_full_blocks = n %/% block_size
    remaining_patients = n %% block_size

  if (remaining_patients %% 2 != 0) {
    stop("Patient number must be even.")
  }

  for (i in 1:num_full_blocks) {
    block = rep(treatments, each = block_size / length(treatments))
    assignments[((i - 1) * block_size + 1):(i * block_size)] = sample(block)
  }

  if (remaining_patients == 2) {
    last_block = sample(c("A", "B"))
    assignments[(n - 1):n] = last_block
  }

  return(assignments)
}


Zhao_generate_data = function(n, m0 = 10, k, beta, isResponseBinary) {
  t = c(permuted_block(2 * m0), rep(NA, n - 2 * m0))
  # Z1=sample(c(1,-1),n,replace = T)
  # Z2=sample(c(1,-1),n,replace = T)
  X = sample(c(1, -1), n, replace = T)
  Z = replicate(k, sample(c(1, -1), n, replace = T))
  # py=1/(1+exp(-(0.5+1*X-0.5*X*t+0.5*Z1+0.5*Z2)))
  if (isResponseBinary) {
  py = 1 / (1 + exp(-(beta[1]*t+ beta[2]*(1-t) + beta[3] * X + beta[4] * X * t +
                        rowSums(beta[5:length(beta)] * Z))))
  py_2m0 = py[1:2 * m0]
  Y = runif(2 * m0)

    Y = c(ifelse(Y < py_2m0, 1, 0), rep(NA, n - 2 * m0))
  } else{
    Y = beta[1]*t+ beta[2]*(1-t) + beta[3] * X + beta[4] * X * t +
      rowSums(beta[5:length(beta)] * Z)+rnorm(n)
  }

  #return(data.frame(cbind(X,Z1,Z2,Y,t)))
  patients=cbind(X, Z, Y, t)
  colnames(patients)[1]="X"
  colnames(patients)[2:(1+k)]=paste0("Z",1:k)
  colnames(patients)[2+k]="Y"
  colnames(patients)[3+k]="t"
  return(patients)
}

d1 = Zhao_generate_data(n=400, beta=c(0.5,-0.5, 0.5, 0.5,0.5,0.5),k = 3, isResponseBinary = F)
nrow(d1)
colnames(d1)
paste0("Z",1:3)
colnames(d1)[1]="X"
colnames(d1)[2:(1+2)]=paste0("Z",1:2)
colnames(d1)[2+2]="Y"
colnames(d1)[3+2]="t"

setequal(d1[1,paste0("Z",1:2)],c(-1,-1))
data_before=d1[1:20,]
lm_before=lm(data_before[,"Y"]~data_before[,"t"]+(1-data_before[,"t"])+
               data_before[,"X"]data_before[,"X"]*data_before[,"t"])
