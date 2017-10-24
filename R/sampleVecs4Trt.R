
## ## ## ## ## ## ## ## ## ## ## ## ## ##
##### -3. Creating vec_dfm_no_alphas #####
## ## ## ## ## ## ## ## ## ## ## ## ## ##



# #' Wrapper function to the creating of the matrix of vectors to evaluate
# #'
# #' @inheritParams argTests
# #' @inheritParams estimateTargets
# #'
# #' @export
sampleVecs4Trt <- function(
  split_data_list,
  max_k_evals
){

  all_trt_info_dfm <- summarizeTrtInfo(split_data_list)

  unique_trt_info_dfm <- getUniqueTrtInfo(
    all_trt_info_dfm
  )
  vec_dfm_no_alphas <- genAllTrtVecs(
    unique_trt_info_dfm,
    max_k_evals
  )
}






# #' Summarize treatment for each cluster
# #'
# #' @inheritParams sampleVecs4Trt
# #'
# #' @export
summarizeTrtInfo <- function(split_data_list){

  group_trts_info <-
    lapply(split_data_list, function(data_list){

      force(data_list)

      treatment_vec <-  data_list$treatment_vec
      len_trt = length(treatment_vec)
      sum_trt = sum(treatment_vec)

      ### defense
      if (len_trt!= nrow(data_list$model_matrix)) {stop('wrong length 1')}
      ### defense

      out <- matrix(c(len_trt, sum_trt), nrow= 1, ncol=2)
      colnames(out) <- c('n','s')

      out
    })

  group_n_s_mat <- as.data.frame(do.call(rbind,group_trts_info),
                                 stringsAsFactors = FALSE)

  ### defense
  if (nrow(group_n_s_mat)!= length(split_data_list)) {stop('wrong length 2')}
  if (ncol(group_n_s_mat)!=2) { stop('wrong length 3')}
  ### defense

  vec_dfm <- unique(group_n_s_mat)

  vec_dfm <- vec_dfm[order(vec_dfm$n,vec_dfm$s),]
  row.names(vec_dfm) <- NULL

  vec_dfm
}


# #' Determine which omegas need to be estimated
# #'
# #' @param all_trt_info_dfm dataframe of number in cluster and sum treated for
# #'   each cluster
# #'
# #' @export
getUniqueTrtInfo <- function(
  all_trt_info_dfm
){
  unique_trt_info_dfm <- by(all_trt_info_dfm, all_trt_info_dfm$n,
                            function(n_data){
                              num_sums <- nrow(n_data)
                              group_size <- n_data$n[[1]]
                              if (num_sums <= group_size){
                                n_data$all_probs <- "none_removed"
                              } else
                                if ( num_sums == (1+group_size)){
                                  n_data$all_probs <- "removed_one"
                                  remove_sum <- floor((num_sums+1)/2)
                                  n_data <- n_data[-remove_sum,]
                                } else { stop('Wrong number of groups maybe?')}
                              n_data
                            })

  unique_trt_info_dfm <- do.call(rbind, unique_trt_info_dfm)
}





# #' Wrapper function to get all treatment vectors to estimate each omega (for one policy)
# #'
# #' @param unique_trt_info_dfm info on all the omegas to estimate
# #' @inheritParams sampleVecs4Trt
# #'
# #' @export
genAllTrtVecs <- function(
  unique_trt_info_dfm,
  max_k_evals
){

  unique_trt_info_dfm$trt_vecs_mats <- NA

  for (ii in 1:nrow(unique_trt_info_dfm)){
    vec_info <- genKTrtVecs_ns(
      n = unique_trt_info_dfm$n[ii],
      s = unique_trt_info_dfm$s[ii],
      max_k = max_k_evals
    )
    vec_mat <- vec_info$all_vecs
    mult_factor <- vec_info$mult_factor
    unique_trt_info_dfm$trt_vecs_mats[ii] <- list(vec_mat)
    unique_trt_info_dfm$mult_factor[ii] <- mult_factor
  }

  unique_trt_info_dfm
}


# #' Generates a certain amount of treatment vectors to evaluate for
# #'
# #' @param  n number of individuals in cluster
# #' @param s sum of number of treated individuals in cluster
# #' @param max_k maximum number of vectors to evaluate
# #'
# #' @export
genKTrtVecs_ns <- function(
  n,
  s,
  max_k
){
  if (max_k==0) {
    max_k <- choose(n,s)
    mult_factor <- 1
  } else
    if (max_k != as.integer(max_k)) {
      stop(paste('max_k needs to be 0 or an integer; argument max_k =',max_k))}

  if (max_k>10 | (max_k==0 & choose(n,s)>10)) {
    message("argument max_k > 10 may result in very slow computation")
  }


  num_samples <- min( choose(n,s),max_k)
  mult_factor <- choose(n,s)/num_samples
  all_vecs <- matrix(NA, ncol = n, nrow = num_samples)

  sample_num = 1
  zero_vec <- rep(0,n)
  while (sample_num <= num_samples) {
    proposal <- sample.int(n,size = s, replace = FALSE )
    proposal_vec <- zero_vec #resetting
    proposal_vec[proposal] <- 1

    if (sample_num > 1 &&
        any(
          apply(all_vecs[1:(sample_num)-1,, drop=FALSE], 1,function(x){
            all.equal(x,proposal_vec, check.attributes = FALSE)==TRUE
          }) ##returns multiple logicals - returns TRUE if duplicated
        ) #returns TRUE if any is duplicated
    ) { ##if any duplicated ==TRUE then skip and start again
      next
    } else {
      all_vecs[sample_num,] <- proposal_vec
      sample_num <- sample_num+1
      next
    }
  }


  out_list <- list(
    all_vecs = all_vecs,
    mult_factor = mult_factor
  )
}
