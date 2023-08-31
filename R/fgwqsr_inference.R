# File containing inference functions for FGWQSR

marginal_prob_outside_direction = function(pollutant_num, mean, sigma, pollutant_direction)
{
  # Arguments:

  # Returns:
  # num_pollutants = ifelse(nrow(sigma) %>% is.na, 1, nrow(sigma))

  marginal_prob = abs(ifelse(pollutant_direction == 1, 0, 1) -
                        stats::pnorm(q = 0,  mean = mean, sd = sigma))

  return(marginal_prob)
}

projection_matrix = function(P) # P is a basis with column vectors corresponding to basis vectors
{
  P%*% solve(t(P) %*% P) %*% t(P)
}

L2_norm = function(x)
{
  return( (x^2) %>% sum %>% sqrt )
}

# function to create list of vectors, where indices denotes the break
# points (inclusive) of the original vector
listify = function(vector, indices)
{
  list = vector(mode = "list", length = length(indices))
  current_index = 1
  for(i in 1:length(indices))
  {
    list[[i]] = vector[current_index: (current_index + indices[i] - 1)]
    current_index = current_index + indices[i]
  }
  return(list)
}



# # function to perform manual projections under H0 and H1
# build_cone_regions = function(effects, cov_mat, group_num, pollutant_num, zero_threshold_cutoff)
# {
#   num_groups = length(effects) # number of groups
#   num_in_each_group = sapply(effects, length)
#
#   effects_vector = unlist(effects) # to have vector of MLE effects
#
#   # group directionality from MLEs
#   pollutant_directions = sapply(effects, sum) %>% sign %>% rep(times = num_in_each_group)
#
#   # generate marginal prob outside group direction under NULL
#   marginal_probs =
#     sapply(1:length(effects_vector),
#            FUN = function(i)
#              return(marginal_prob_outside_direction(pollutant_num = i,
#                                                     mean = effects_vector[i],
#                                                     sigma = sqrt(cov_mat[i,i]),
#                                                     pollutant_direction = pollutant_directions[i])))
#
#   # generate cone regions -- group of interest should be zero'd out.
#
#   region_vector = rep("R1", sum(num_in_each_group)) # initialize to R1
#
#   # find which pollutants on boundary by zero_threshold_cutoff
#   boundary_cases = which(marginal_probs >= zero_threshold_cutoff)
#
#   # change boundary cases to half line cone for pollutants on a boundary
#   for(i in boundary_cases)
#   {
#     region_vector[i] = ifelse(pollutant_directions[i] == -1,
#                               "(-Inf,0]", "[0,Inf)")
#   }
#
#   # if there is a group with one pollutant (ie, in the case of H0, set to R1)
#
#   one_pollutant_groups = which(num_in_each_group == 1)
#
#   if(length(one_pollutant_groups) > 0) # if there is at least one 1 pollutant group
#   {
#     for(i in one_pollutant_groups)
#     {
#       pollutant_loc = ifelse(i == 1, 1, num_in_each_group[1:i-1] %>% sum(1))
#       # change to R1 cone
#       region_vector[pollutant_loc] = "R1"
#
#       if(pollutant_loc %in% boundary_cases)
#       {
#         boundary_cases = boundary_cases[-which(boundary_cases == pollutant_loc)]
#       }
#     }
#   }
#
#   # boolean vector if each pollutant is on boundary or not
#   boundary_cases_boolean = 1:sum(num_in_each_group) %in% boundary_cases
#
#   # change cone to sign constrained region for nuisance pollutants with groups that have params all on boundary.
#   for(i in (1:num_groups))
#   {
#     # pollutant numbers for each group
#     pollutant_locs_by_group = ifelse(i == 1, 1, num_in_each_group[1:i-1] %>% sum(1)) :
#       sum(num_in_each_group[1:i])
#     # get subvector for group of whether pollutants on boundary
#     group_subvector_boundary = boundary_cases_boolean[pollutant_locs_by_group]
#
#     # if number of pollutants on boundary is same as number of pollutants in the group, set SC cone
#     if(sum(group_subvector_boundary) == num_in_each_group[i])
#     {
#       region_vector[pollutant_locs_by_group] = "SC"
#     }
#   }
#
#   return(region_vector)
# }
#
# generate_group_LRT = function(Z, effects, cov_mat,
#                               affected_pollutants,
#                               cone_regions,
#                               cone_regions_list,
#                               nuis_groups_cat,
#                               group_directions,
#                               num_in_each_group,
#                               group_num)
# {
#   # number pollutants full mod
#   p = effects %>% unlist %>% length
#
#   # basis for R^p -- columns will serve as basis vectors in original space
#   full_basis = diag(p)
#
#   # initialize for later
#   remove_basis_vectors_nuisance = numeric()
#
#   # find basis under H0 -- need to iterate over SC and half plane nuisance groups.
#   # interior R^p_k groups will always have their basis vectors included
#
#   # first, iterate over half plane groups
#
#   if(length(nuis_groups_cat$half_plane_nuis_groups) > 0) #if half plane nuisance groups
#   {
#     HLP_indecies = which(cone_regions %in% c("[0,Inf)", "(-Inf,0]"))
#     HLP_group_directions = rep(group_directions, times = num_in_each_group)[HLP_indecies]
#
#     HLP_direction_violators = HLP_indecies[which(sign(Z[HLP_indecies]) != HLP_group_directions)]
#
#     remove_basis_vectors_nuisance = HLP_direction_violators
#   }
#
#   # now consider nuisance groups that have SC cones
#   # if we have nuisance groups with SC cone
#
#   basis_proposal_list = list() # stores candidate projections for groups with SC
#
#   # if there are nuisance groups with SC cones
#   # group of interest will not be included in sc_nuis_groups
#   if(length(nuis_groups_cat$sc_nuis_groups) > 0)
#   {
#     basis_sc_list_H0 = list()
#     basis_sc_list_H0_index = 1
#
#     projected_SC = F # stores if we project a SC group
#
#     for(i in nuis_groups_cat$sc_nuis_groups) # iterate over SC groups
#     {
#       sc_poll_indices = cone_regions_list[[i]]$pollutant_indices
#       pollutant_signs = Z[sc_poll_indices] %>% sign
#
#       # if not all signs in the group are the same
#       if(abs(pollutant_signs %>% sum) != length(sc_poll_indices))
#       {
#         # remove basis vectors from candidates
#         remove_basis_vectors_nuisance = append(remove_basis_vectors_nuisance, sc_poll_indices)
#
#         # find locations of negative and positive pollutants
#         pos_polls = which(pollutant_signs == 1)
#
#         neg_polls = which(pollutant_signs == -1)
#
#         # append to list the two possible basis vector configurations
#         # either to include basis vectors corresponding to positive
#         # pollutants or to negative pollutants
#
#         basis_sc_list_H0[[basis_sc_list_H0_index]] =
#           list(pos = sc_poll_indices[pos_polls],
#                neg = sc_poll_indices[neg_polls])
#
#         basis_sc_list_H0_index = basis_sc_list_H0_index +1
#         projected_SC = T # this is so we save that we need to project
#       }
#     }
#
#     if(projected_SC == T) # if we projected a SC group
#     {
#       # find permutations of basis vector inclusions
#
#       region_permutations = permutations(n = 2, r = length(basis_sc_list_H0) ,
#                                          v = c("neg","pos"), repeats.allowed = T)
#
#       for(i in 1: nrow(region_permutations))
#       {
#         basis = numeric()
#         for(j in 1: ncol(region_permutations))
#         {
#           basis_indices = basis_sc_list_H0[[j]][[region_permutations[i,j]]]
#           basis = append(basis, basis_indices)
#         }
#         basis_proposal_list[[i]] = basis
#       }
#     }
#   }
#
#   # add affected pollutants to remove from basis under H0
#   remove_basis_vectors_nuisance = append(remove_basis_vectors_nuisance,
#                                          affected_pollutants)
#
#   # now begin formulating projections
#   eigen_info = eigen(cov_mat, symmetric = T)
#   P = eigen_info$vectors
#   lambda = eigen_info$values
#   # transformation to N(lambda^-.5P^T * mu, I) space
#   trans_mat = diag(lambda^(-.5)) %*% t(P)
#   Z_tilde = matrix(tcrossprod(trans_mat, matrix(Z, nrow = 1)), nrow = 1)
#
#   trans_projection_H0 = 0
#
#   if(length(basis_proposal_list) == 0) # If we did not perform any SC projections
#   {
#     # quick way to perform the matrix multiplication of trans_mat %*% basis
#     exclude_indices = remove_basis_vectors_nuisance
#     basis_transformed_H0 = trans_mat[, -exclude_indices]
#     proj_mat_H0 = projection_matrix(basis_transformed_H0)
#     trans_projection_H0 = tcrossprod(proj_mat_H0, Z_tilde)
#
#   }else # there were SC projected pollutants
#   {
#     # iterate over proposal bases
#     proposal_projection_list_H0 = vector(mode = "list",
#                                       length(basis_proposal_list))
#
#     for(i in 1: length(basis_proposal_list))
#     {
#       # denote indices to exclude
#       exclude_indices =  setdiff(remove_basis_vectors_nuisance,basis_proposal_list[[i]])
#       # create basis
#       basis_transformed_H0 = trans_mat[, -exclude_indices]
#
#       # create projection matrix
#       proj_mat_H0 = projection_matrix(basis_transformed_H0)
#
#       # perform candidate projection
#       projection_H0 = tcrossprod(proj_mat_H0, Z_tilde)
#
#       # save projection and norm
#       proposal_projection_list_H0[[i]] = list(projection = projection_H0,
#                                            norm = L2_norm(projection_H0))
#     }
#
#     max_proj = which.max(sapply(proposal_projection_list_H0, "[[", 2))
#
#     trans_projection_H0 = proposal_projection_list_H0[[max_proj]]$projection
#   }
#
#   # find projection under H_1
#
#   trans_projection_H1 = 0 #initialize
#   # indices for group of interest
#   sc_poll_indices = affected_pollutants
#
#   # signs of affected pollutants
#   pollutant_signs = Z[sc_poll_indices] %>% sign
#
#   # find locations of negative and positive pollutants
#   pos_polls = affected_pollutants[which(pollutant_signs == 1)]
#
#   neg_polls = affected_pollutants[which(pollutant_signs == -1)]
#
#
#   # create new proposal bases from previous proposals
#   basis_proposal_list_H1_1 = list()
#   basis_proposal_list_H1_2 = list()
#   basis_proposal_list_H1 = list()
#
#   if(length(basis_proposal_list)>0) # if there are sign constrained nuisance groups
#   {
#     if(length(pos_polls) > 0) # if there are positive pollutants in group of interest
#     {
#       basis_proposal_list_H1_1 = lapply(basis_proposal_list,
#                                         FUN = function(x) return(c(x,pos_polls )))
#     }
#     if(length(neg_polls) > 0) # if there are negative pollutants in group of interest
#     {
#       basis_proposal_list_H1_2 = lapply(basis_proposal_list,
#                                         FUN = function(x) return(c(x,neg_polls )))
#     }
#     basis_proposal_list_H1 = append(basis_proposal_list_H1_1, basis_proposal_list_H1_2)
#   }else # if there are no sign constrained nuisance groups
#   {
#     if(length(pos_polls) > 0) # if there are positive pollutants in group of interest
#     {
#       basis_proposal_list_H1 = append(basis_proposal_list_H1, list(pos = pos_polls))
#     }
#     if(length(neg_polls) > 0) # if there are negative pollutants in group of interest
#     {
#       basis_proposal_list_H1 = append(basis_proposal_list_H1, list(neg = neg_polls) )
#     }
#
#   }
#
#
#   proposal_projection_list_H1 = vector(mode = "list",
#                                        length(basis_proposal_list_H1))
#
#   for(i in 1: length(basis_proposal_list_H1))
#   {
#     # denote indices to exclude
#     exclude_indices =  setdiff(remove_basis_vectors_nuisance,
#                                basis_proposal_list_H1[[i]])
#     basis_transformed_H1 = NULL
#
#     if(length(exclude_indices) != 0) # if there are excluded indices
#     {
#     # create basis
#     basis_transformed_H1 = trans_mat[, -exclude_indices]
#     }else # if there are no excluded indices
#     {
#       basis_transformed_H1 = trans_mat
#     }
#
#     # create projection matrix
#     proj_mat_H1 = projection_matrix(basis_transformed_H1)
#
#     # perform candidate projection
#     projection_H1 = tcrossprod(proj_mat_H1, Z_tilde)
#
#     # save projection and norm
#     proposal_projection_list_H1[[i]] = list(projection = projection_H1,
#                                             norm = L2_norm(projection_H1))
#   }
#
#   max_proj = which.max(sapply(proposal_projection_list_H1, "[[", 2))
#
#   trans_projection_H1 = proposal_projection_list_H1[[max_proj]]$projection
#
#
#
#   # compute LRT -- mvn log likelihood does not use -1/2.
#
#   LRT = mvn_log_likelihood(mean = c(trans_projection_H0), sample = c(Z_tilde),
#                            inv_cov_mat = diag(p)) -
#     mvn_log_likelihood(mean = c(trans_projection_H1), sample = c(Z_tilde),
#                        inv_cov_mat = diag(p))
#
#   return(LRT)
# }
#
# generate_SPLRT = function(Z, effects, cov_mat,
#                           affected_pollutants,
#                           cone_regions,
#                           cone_regions_list,
#                           nuis_groups_cat,
#                           group_directions,
#                           num_in_each_group,
#                           group_num,
#                           pollutant_num)
# {
#   # number pollutants full mod
#   p = effects %>% unlist %>% length
#
#   # basis for R^p -- columns will serve as basis vectors in original space
#   full_basis = diag(p)
#
#   # initialize for later
#   remove_basis_vectors_nuisance = numeric(p)
#
#   # find basis under H0 -- need to iterate over SC and half plane nuisance groups.
#   # interior R^p_k groups will always have their basis vectors included
#
#   # first, iterate over half plane groups
#
#   if(length(nuis_groups_cat$half_plane_nuis_groups) > 0) #if half plane nuisance groups
#   {
#     HLP_indecies = which(cone_regions %in% c("[0,Inf)", "(-Inf,0]"))
#     HLP_group_directions = rep(group_directions, times = num_in_each_group)[HLP_indecies]
#
#     HLP_direction_violators = HLP_indecies[which(sign(Z[HLP_indecies]) != HLP_group_directions)]
#
#     remove_basis_vectors_nuisance = HLP_direction_violators
#   }
#
#   # now consider nuisance groups that have SC cones
#   # if we have nuisance groups with SC cone
#
#   basis_proposal_list = list() # stores candidate projections for groups with SC
#   projected_SC = F # stores if we project a SC group
#   basis_sc_list_H0 = list()
#   basis_sc_list_H0_index = 1
#
#   # if there are nuisance groups with SC cones
#   # group of interest will not be included in sc_nuis_groups
#   if(length(nuis_groups_cat$sc_nuis_groups) > 0)
#   {
#
#     for(i in nuis_groups_cat$sc_nuis_groups) # iterate over SC groups
#     {
#       sc_poll_indices = cone_regions_list[[i]]$pollutant_indices
#       pollutant_signs = Z[sc_poll_indices] %>% sign
#
#       # if not all signs in the group are the same
#       if(abs(pollutant_signs %>% sum) != length(sc_poll_indices))
#       {
#         # remove basis vectors from candidates
#         remove_basis_vectors_nuisance = append(remove_basis_vectors_nuisance, sc_poll_indices)
#
#         # find locations of negative and positive pollutants
#         pos_polls = which(pollutant_signs == 1)
#
#         neg_polls = which(pollutant_signs == -1)
#
#         # append to list the two possible basis vector configurations
#         # either to include basis vectors corresponding to positive
#         # pollutants or to negative pollutants
#
#         basis_sc_list_H0[[basis_sc_list_H0_index]] =
#           list(pos = sc_poll_indices[pos_polls],
#                neg = sc_poll_indices[neg_polls])
#
#         basis_sc_list_H0_index = basis_sc_list_H0_index +1
#         projected_SC = T # this is so we save that we need to project
#       }
#     }
#   }
#
#   # if the group from where the pollutant of interest resides has SC cone
#   if(cone_regions_list[[group_num]]$regions[1] == "SC")
#   {
#     sc_poll_indices = cone_regions_list[[group_num]]$pollutant_indices[-pollutant_num]
#     pollutant_signs = Z[sc_poll_indices] %>% sign
#
#     # if not all signs in the group are the same
#     if(abs(pollutant_signs %>% sum) != length(sc_poll_indices))
#     {
#       # remove basis vectors from candidates
#       remove_basis_vectors_nuisance = append(remove_basis_vectors_nuisance, sc_poll_indices)
#
#       # find locations of negative and positive pollutants
#       pos_polls = which(pollutant_signs == 1)
#
#       neg_polls = which(pollutant_signs == -1)
#
#       # append to list the two possible basis vector configurations
#       # either to include basis vectors corresponding to positive
#       # pollutants or to negative pollutants
#
#       basis_sc_list_H0[[basis_sc_list_H0_index]] = list(pos = sc_poll_indices[pos_polls],
#                                                         neg = sc_poll_indices[neg_polls])
#       projected_SC = T # this is so we save that we need to project
#     }
#
#   }
#
#   if(projected_SC == T) # if we projected a SC group
#   {
#     # find permutations of basis vector inclusions
#
#     region_permutations = permutations(n = 2, r = length(basis_sc_list_H0) ,
#                                        v = c("neg","pos"), repeats.allowed = T)
#
#     for(i in 1: nrow(region_permutations))
#     {
#       basis = numeric()
#       for(j in 1: ncol(region_permutations))
#       {
#         basis_indices = basis_sc_list_H0[[j]][[region_permutations[i,j]]]
#         basis = append(basis, basis_indices)
#       }
#       basis_proposal_list[[i]] = basis
#     }
#   }
#
#   # add affected pollutant to remove from basis under H0
#   remove_basis_vectors_nuisance = append(remove_basis_vectors_nuisance,
#                                          affected_pollutants)
#
#   # now begin formulating projections
#   eigen_info = eigen(cov_mat, symmetric = T)
#   P = eigen_info$vectors
#   lambda = eigen_info$values
#   # transformation to N(lambda^-.5P^T * mu, I) space
#   trans_mat = diag(lambda^(-.5)) %*% t(P)
#   Z_tilde = matrix(tcrossprod(trans_mat, matrix(Z, nrow = 1)), nrow = 1)
#
#   basis_H0 = numeric()
#   trans_projection_H0 = 0
#
#   if(length(basis_proposal_list) == 0) # If we did not perform any SC projections
#   {
#     # quick way to perform the matrix multiplication of trans_mat %*% basis
#     exclude_indices =  remove_basis_vectors_nuisance
#     basis_transformed_H0 = trans_mat[, -exclude_indices]
#     proj_mat_H0 = projection_matrix(basis_transformed_H0)
#
#     trans_projection_H0 = tcrossprod(proj_mat_H0, Z_tilde)
#     basis_H0 = setdiff(1:p, remove_basis_vectors_nuisance)
#
#   }else # there were SC projected pollutants
#   {
#     # iterate over proposal bases
#     proposal_projection_list = vector(mode = "list",
#                                       length(basis_proposal_list))
#
#     for(i in 1: length(basis_proposal_list))
#     {
#       # denote indices to exclude
#       exclude_indices = c(affected_pollutants, setdiff(remove_basis_vectors_nuisance,basis_proposal_list[[i]]))
#       # create basis
#       basis_transformed_H0 = trans_mat[, -exclude_indices]
#
#       # create projection matrix
#       proj_mat_H0 = projection_matrix(basis_transformed_H0)
#
#       # perform candidate projection
#       projection_H0 = tcrossprod(proj_mat_H0, Z_tilde )
#
#       # save projection and norm
#       proposal_projection_list[[i]] = list(projection = projection_H0,
#                                            norm = L2_norm(projection_H0),
#                                            basis = setdiff(1:p, exclude_indices))
#     }
#
#     max_proj = which.max(sapply(proposal_projection_list, "[[", 2))
#
#     trans_projection_H0 = proposal_projection_list[[max_proj]]$projection
#     basis_H0 = proposal_projection_list[[max_proj]]$basis
#   }
#
#   # find projection under H_1
#
#   trans_projection_H1 = 0 #initialize
#
#   # if the group where the pollutant of interest resides is NOT SC
#   if(cone_regions_list[[group_num]]$regions[1] != "SC")
#   {
#     # if the pollutant has Z in the same direction as the group
#     if(sign(Z[affected_pollutants]) == group_directions[group_num])
#     {
#       basis_H1 = append(basis_H0,affected_pollutants )
#
#       basis_transformed_H1 = trans_mat[, basis_H1]
#
#       # create projection matrix
#       proj_mat_H1 = projection_matrix(basis_transformed_H1)
#
#       # perform candidate projection
#       trans_projection_H1 = tcrossprod(proj_mat_H1, Z_tilde )
#     }else # if the pollutant has Z NOT in the same direction as the group
#     {
#       trans_projection_H1 = trans_projection_H0 # same as under H0
#     }
#   }else # if the group where pollutant of interest resides has SC cone
#   {
#
#     sc_poll_indices = cone_regions_list[[group_num]]$pollutant_indices
#     pollutant_signs = Z[sc_poll_indices] %>% sign
#
#     basis_H0 = setdiff(basis_H0,sc_poll_indices)
#
#     # if not all signs in the group are the same
#     if(abs(pollutant_signs %>% sum) != length(sc_poll_indices))
#     {
#
#       # find locations of negative and positive pollutants
#       pos_polls = which(pollutant_signs == 1)
#
#       neg_polls = which(pollutant_signs == -1)
#
#       # create two candidate projections
#       basis_transformed_H1_1 = trans_mat[, union(basis_H0, sc_poll_indices[pos_polls])]
#       proj_mat_H1_1 = projection_matrix(basis_transformed_H1_1)
#       trans_projection_H1_1 = tcrossprod(proj_mat_H1_1, Z_tilde)
#
#       basis_transformed_H1_2 = trans_mat[, union(basis_H0, sc_poll_indices[neg_polls])]
#       proj_mat_H1_2 = projection_matrix(basis_transformed_H1_2)
#       trans_projection_H1_2 = tcrossprod(proj_mat_H1_2, Z_tilde)
#
#       trans_projection_H1_list = list(trans_projection_H1_1, trans_projection_H1_2)
#
#       trans_projection_H1 = trans_projection_H1_list[[which.max(sapply(trans_projection_H1_list, L2_norm))]]
#     }
#   }
#
#   # compute LRT -- mvn log likelihood does not use -1/2.
#
#   LRT = mvn_log_likelihood(mean = c(trans_projection_H0), sample = c(Z_tilde),
#                            inv_cov_mat = diag(p)) -
#     mvn_log_likelihood(mean = c(trans_projection_H1), sample = c(Z_tilde),
#                        inv_cov_mat = diag(p))
#
#   return(LRT)
# }
#
# mvn_lrt_simulation= function(effects, cov_mat,
#                              group_num, pollutant_num = NA,
#                              num_sims = 10000, zero_threshold_cutoff = .5,
#                              seed = 2023,
#                              cores = availableCores())
# {
#   # number of groups and number of elements in each group
#   num_groups = length(effects)
#   num_in_each_group = sapply(effects, length)
#
#   # simulate the MVN 1 sample realizations
#   set.seed(seed)
#   mvn_samples = rmvnorm(num_sims, mean = effects %>% unlist, sigma = cov_mat)
#
#   lrts = numeric(num_sims) # initialize
#
#   # generate projection under H0 and H1
#   # pollutants of interest should have effects set to 0!
#   cone_regions = build_cone_regions(effects = effects, cov_mat = cov_mat,
#                                     group_num = group_num,
#                                     pollutant_num = pollutant_num,
#                                     zero_threshold_cutoff = zero_threshold_cutoff)
#
#   # make in list by groups to traverse easier
#   cone_regions_list = listify(cone_regions, num_in_each_group)
#
#   # find affected pollutants
#   affected_pollutants = ifelse(group_num == 1, 1, num_in_each_group[1:(group_num-1)] %>% sum(1))
#   if(is.na(pollutant_num)) # if it is a group LRT
#   {
#     # find pollutant group to set to 0 under null
#     affected_pollutants = affected_pollutants:
#       (affected_pollutants + num_in_each_group[group_num] - 1)
#   }else # if it is a single pollutant LRT
#   {
#     # find pollutant to set to 0 under null
#     affected_pollutants = affected_pollutants + pollutant_num - 1
#   }
#
#   # Nuisance group identification
#   # denote nuisance groups that have all effects on interior-- all R^1 cone
#   interior_nuis_groups = sapply(cone_regions_list, FUN = function(x)
#   {
#     times = length(x)
#     return( sum(x == rep("R1", times)) == times )
#   }) %>% which
#
#   # denote nuisance groups with some effects on interior
#   # just need to check if has [0,Inf) or (-Inf, 0], since by definition
#   # cant all be halflines in group -- otherwise would be SC cone
#   half_plane_nuis_groups = sapply(cone_regions_list, FUN = function(x)
#   {
#     return( sum(grepl("[0,Inf)", x, fixed = T) | grepl("(-Inf,0]", x, fixed = T)) >0)
#   }) %>% which
#
#   # use set diff to find SC nuisance groups -- remove group of interest
#   sc_nuis_groups = setdiff((1:num_groups)[-group_num],
#                            c(interior_nuis_groups, half_plane_nuis_groups))
#
#   nuis_groups_cat = list(interior_nuis_groups = interior_nuis_groups,
#                          half_plane_nuis_groups = half_plane_nuis_groups,
#                          sc_nuis_groups = sc_nuis_groups)
#
#   cone_regions_list = listify_cone_regions(cone_regions,num_in_each_group, half_plane_nuis_groups)
#
#   group_directions = sapply(effects, FUN = function(x) { return(sum(x) %>% sign)})
#
#   # perform projections and compute LRT
#
#   # if we are performing a 1 group LRT
#   if(is.na(pollutant_num))
#   {
#     plan(multisession, workers = cores)
#
#     lrts = future_apply(mvn_samples, MARGIN = 1,
#                         FUN = generate_group_LRT,
#                         effects = effects,
#                         cov_mat = cov_mat,
#                         affected_pollutants = affected_pollutants,
#                         cone_regions = cone_regions,
#                         cone_regions_list = cone_regions_list,
#                         nuis_groups_cat = nuis_groups_cat,
#                         group_directions = group_directions,
#                         num_in_each_group = num_in_each_group,
#                         group_num = group_num)
#
#     plan(sequential)
#   }else
#   {
#     plan(multisession, workers = cores)
#
#     lrts = future_apply(mvn_samples, MARGIN = 1,
#                         FUN = generate_SPLRT,
#                         effects = effects,
#                         cov_mat = cov_mat,
#                         affected_pollutants = affected_pollutants,
#                         cone_regions = cone_regions,
#                         cone_regions_list = cone_regions_list,
#                         nuis_groups_cat = nuis_groups_cat,
#                         group_directions = group_directions,
#                         num_in_each_group = num_in_each_group,
#                         group_num = group_num,
#                         pollutant_num = pollutant_num)
#
#     plan(sequential)
#   }
#
#   return(lrts)
# }
#
#
# mvn_log_likelihood = function(mean, sample, inv_cov_mat)
# {
#   return(crossprod(matrix(sample - mean, ncol = 1), inv_cov_mat) %>%
#            tcrossprod(matrix(sample - mean, ncol = 1) %>% t))
# }
#
# # function to create list of vectors, where indices denotes the break
# # points (inclusive) of the original vector
# listify = function(vector, indices)
# {
#   list = vector(mode = "list", length = length(indices))
#   current_index = 1
#   for(i in 1:length(indices))
#   {
#     list[[i]] = vector[current_index: (current_index + indices[i] - 1)]
#     current_index = current_index + indices[i]
#   }
#   return(list)
# }
#
# listify_cone_regions = function(cone_regions, indices, half_plane_nuis_groups)
# {
#   list = vector(mode = "list", length = length(indices))
#   current_index = 1
#   for(i in 1:length(indices))
#   {
#     list[[i]]$regions = cone_regions[current_index: (current_index + indices[i] - 1)]
#     list[[i]]$pollutant_indices = current_index: (current_index + indices[i] - 1)
#     current_index = current_index + indices[i]
#   }
#   return(list)
# }
#

################################################################################
# MVN Multiple Optim Code for Inference


# for single pollutant LRT, modify function

generate_optim_limits = function(effects, cov_mat)
{
  # these are not passed to the function because we change effects and cov_mat depending
  # on H0 and H1, this is the easier implementation to run the function for H0 and H1
  # seperatley

  num_groups = effects %>% length
  num_in_each_group = sapply(effects, length)

  # group directionality from MLEs
  pollutant_directions = sapply(effects, sum) %>% sign %>% rep(times = num_in_each_group)

  # Generate marginal probs that a pollutant is outside of its group direction, will determine whether
  # it is close to the boundary of the constrained region
  # first col has pollutant num, second col has prob outside direction

  #if cov_mat is a single number
  variances = diag(cov_mat)

  marginal_probs = (1:sum(num_in_each_group)) %>% sapply(FUN = function(x)
  {
    return(c(x,
             marginal_prob_outside_direction(x, mean = (effects %>% unlist)[x],
                                             sigma = sqrt(variances[x]),
                                             pollutant_direction = pollutant_directions[x])))
  }) %>% t


  # set initial cones to R1 cones -- we will modify these
  lower = rep(-Inf, length(unlist(effects)))
  upper = rep(Inf, length(unlist(effects)))

  # determine which pollutants are on the boundary by zero_threshold_cutoff
  # we already evaluate the ZTC cutoff in the mvn_simulation function
  boundary_cases = which(unlist(effects) == 0)

  # apply the boundary cone
  for(i in boundary_cases)
  {
    lower[i] = ifelse(pollutant_directions[marginal_probs[i,1]] == -1, -Inf, 0)
    upper[i] = ifelse(pollutant_directions[marginal_probs[i,1]] == -1, 0, Inf)
  }

  # Change cone to sign constraint if all params in a group are on the boundary

  # Boolean vector if each pollutant is near(on) boundary or not on boundary by zero_threshold_cutoff criteria
  boundary_cases_boolean = 1:sum(num_in_each_group) %in% marginal_probs[boundary_cases,1]

  groups_on_boundary = numeric(length = num_groups) # initialize for loop -- 0 means not all pollutants in group on boundary
  # 1 means all pollutants in group are on boundary

  for(i in (1:num_groups))  # iterate over groups, will check if each group has all pollutants on boundary
  {
    # get subvector for group of whether pollutants on boundary
    group_subvector_boundary = boundary_cases_boolean[ ifelse(i == 1, 1, num_in_each_group[1:i-1] %>% sum(1)) :
                                                         sum(num_in_each_group[1:i])]

    # if sum of all boolean values for the group equal the group size, then all pollutants in the
    # group are on boundary -- make them have sign constrained cone
    groups_on_boundary[i] = ifelse(sum(group_subvector_boundary) == num_in_each_group[i], 1, 0)
  }

  # find which groups need sign constrained cone
  groups_on_boundary = (groups_on_boundary == 1) %>% which

  # if a group contains one pollutant, we need to ensure it has an R1 cone
  # and potentially remove it from the list of groups on the boundary

  one_pollutant_groups = which(num_in_each_group == 1) # which groups have 1 element

  for(i in 1:length(one_pollutant_groups))
  {
    # one pollutant group to modfiy limits
    one_poll_group = one_pollutant_groups[i]

    # find index number in limit vector -- note it is just one element
    pollutant_loc = ifelse(one_poll_group == 1,1,
                           num_in_each_group[1:(one_poll_group-1)] %>% sum(1))
    # set R1 cone
    lower[pollutant_loc] = -Inf
    upper[pollutant_loc] = Inf

  }

  # remove one pollutant groups from vector of groups on the boundary

  groups_on_boundary = base::setdiff(groups_on_boundary,one_pollutant_groups)

  # now we find the different areas of optimization

  optimization_limits_list = list()

  if(length(groups_on_boundary) == 0) # if no nuisance groups on boundary
  {
    optimization_limits_list = list(list(lower = lower,
                                         upper = upper))

  }else # if there are nuisance groups with entire group on boundary
  {
    # create permutations of sign constrained regions for groups with sign constrainted approx cone
    permutations_regions = gtools::permutations(n = 2, r = groups_on_boundary %>% length,
                                        v = c(-1,1), repeats.allowed = T)

    # add sign constrained cones to limits of optimization
    for(i in 1:nrow(permutations_regions)) # iterate over permutation region optimization number
    {
      lower_permute = lower # initialize these, will be changed
      upper_permute = upper

      for(j in 1: ncol(permutations_regions)) # iterate over groups that need sign constrained region
      {
        # group number to apply sign constraint
        boundary_group = groups_on_boundary[j]

        # find location of affected pollutants
        # this line is complicated because finding the group before -- need to account for discontinuity from group removal of vector
        pollutant_locs = ifelse(boundary_group == 1,1,
                                num_in_each_group[1:(boundary_group-1)] %>% sum(1)) :
          sum(num_in_each_group[1:boundary_group])

        # apply sign constraint  -- -1 from permutation indicates (-Inf, 0], 1 indicates [0,Inf)
        if(permutations_regions[i,j] == -1)
        {
          lower_permute[pollutant_locs] = -Inf
          upper_permute[pollutant_locs] = 0
        }else
        {
          lower_permute[pollutant_locs] = 0
          upper_permute[pollutant_locs] = Inf
        }
      }
      optimization_limits_list[[i]] = list(lower = lower_permute,
                                           upper = upper_permute) # add the new optimization limits to the list of optimizations
    }
  }

  return(optimization_limits_list)
}


# generate intitial values that are inside constrained optimization regions
generate_initial_vals = function(optimization_lower, optimization_upper)
{
  initial_vals = numeric(optimization_lower %>% length)

  for(i in 1: length(optimization_lower))
  {
    if(optimization_lower[i] == 0) # if the lower bound is 0, then it is from [0, Inf)
    {
      initial_vals[i] = 1

    }else if(optimization_upper[i] == 0) # if upper bound is 0, then it is from (-Inf, 0]
    {
      initial_vals[i] = -1
    }else # if the variable is in (-Inf, Inf)
    {
      initial_vals[i] = 0
    }

  }
  return(initial_vals)
}


mvn_lrt_simulation= function(effects, cov_mat, group_num, pollutant_num = NA,
                             num_sims = 10000, zero_threshold_cutoff = .5,
                             cores = future::availableCores(), seed = 2023)
{
  num_groups = length(effects) # number of groups
  num_in_each_group = sapply(effects, length) # number of pollutants in each group

  # Generate marginal probs that a pollutant is outside of its group direction, will determine whether
  # it is close to the boundary of the constrained region
  # first col has pollutant num, second col has prob outside direction

  variances = diag(cov_mat)
  # group directional from MLEs
  pollutant_directions = sapply(effects, sum) %>% sign %>% rep(times = num_in_each_group)

  marginal_probs = (1:sum(num_in_each_group)) %>% sapply(FUN = function(x)
  {
    return(
      marginal_prob_outside_direction(x, mean = (effects %>% unlist)[x],
                                      sigma = sqrt(variances[x]),
                                      pollutant_direction = pollutant_directions[x]))
  })

  # set effects of parameters with marginal prob < ZTC to 0
  effects = ifelse(marginal_probs >= as.numeric(zero_threshold_cutoff),
                          0,unlist(effects))
  effects  = listify(effects , num_in_each_group)

  # simulate the MVN 1 sample realizations
  set.seed(seed) # this is so if we repeat inference, we should get same pvalue
  mvn_samples = mvtnorm::rmvnorm(num_sims, mean = effects %>% unlist, sigma = cov_mat)
  inv_cov_mat = pracma::pinv(cov_mat) # used in optim, use moore pentrose psuedoinverse

  lrts = numeric(num_sims) # initialize

  # if there is only one group and we are performing group LRT, we find likelihood
  # when all params set to 0
  if((num_groups == 1) && (is.na(pollutant_num)))
  {
    # initialize the limits of optimization under H1
    optimization_limits_list_H1 = list(
      list(lower = rep(-Inf,num_in_each_group ), upper = rep(0, num_in_each_group)),
      list(lower = rep(0,num_in_each_group ), upper = rep(Inf, num_in_each_group)))

    # perform optimization under H_0 and H_1

    future::plan(future::multisession, workers = cores)

    ll_H0 = -1*future.apply::future_apply(mvn_samples, MARGIN = 1, FUN = function(x)
    {
      return(mvn_log_likelihood(mean = rep(0, num_in_each_group), sample = x, inv_cov_mat = inv_cov_mat))
    })

    ll_H1 = future.apply::future_apply(mvn_samples, MARGIN = 1,
                         FUN = lrt_optimization_H1,
                         inv_cov_mat = inv_cov_mat,
                         optimization_limits_list = optimization_limits_list_H1)

    future::plan(future::sequential)

    lrts = (ll_H1 - ll_H0) # our ll did not divide by the 2, so we leave it out

  }else # if we are performing group LRT with more than one group OR single pollutant LRT
  {

    # we modify this in the if else
    affected_pollutants = ifelse(group_num == 1, 1, num_in_each_group[1:(group_num-1)] %>% sum(1))

    if(is.na(pollutant_num)) # if it is a group LRT
    {
      # find pollutant group to set to 0 under null
      affected_pollutants = affected_pollutants:(affected_pollutants + num_in_each_group[group_num] - 1)
    }else # if it is a single pollutant LRT
    {
      # find pollutant to set to 0 under null
      affected_pollutants = affected_pollutants + pollutant_num - 1
    }

    # remove pollutants set to 0 from cov_mat under null -- this will make it so that we
    # do not evaluate the marginal probabilities outside interior cone for these pollutants
    cov_mat_H0 = cov_mat[-affected_pollutants, -affected_pollutants]

    effects_H0 = effects # we modify this below

    if(is.na(pollutant_num)) # if it is a group lrt
    {
      # remove vector of effects for the group from effects vector under H0
      effects_H0 = effects_H0[-group_num]
    }else # if it is a single pollutant LRT
    {
      # remove pollutant from effects vector under H0
      effects_H0[[group_num]] = effects_H0[[group_num]][-pollutant_num]
    }

    optimization_limits_list_H0 = generate_optim_limits(effects = effects_H0,
                                                        cov_mat = cov_mat_H0)

    optimization_limits_list_H1 = generate_optim_limits(effects = effects,
                                                        cov_mat = cov_mat)

    # perform optimization

    future::plan(future::multisession, workers = cores)

    ll_H0 =  future.apply::future_apply(mvn_samples, MARGIN = 1,
                         FUN = lrt_optimization_H0,
                         inv_cov_mat = inv_cov_mat,
                         optimization_limits_list = optimization_limits_list_H0,
                         affected_pollutants = affected_pollutants)

    ll_H1 = future.apply::future_apply(mvn_samples, MARGIN = 1,
                         FUN = lrt_optimization_H1,
                         inv_cov_mat = inv_cov_mat,
                         optimization_limits_list = optimization_limits_list_H1)

    future::plan(future::sequential)

    lrts = (ll_H1 - ll_H0) # our likelihood function did not divide by 2, so we exclude the 2

  }

  return(lrts)
}

# lrt optimizations H0 and H1 call different mvn log likelihood functions
# (H0 appends 0s for affected pollutants in mvn function)
lrt_optimization_H0 = function(mvn_sample, inv_cov_mat,optimization_limits_list, affected_pollutants)
{
  ll_vals = numeric(length(optimization_limits_list)) # stores the log likelihood

  for(i in 1: length(optimization_limits_list))
  {

    initial_values = generate_initial_vals(optimization_lower = optimization_limits_list[[i]]$lower,
                                           optimization_upper = optimization_limits_list[[i]]$upper)
    # this try catch is because of the following error:
    # Error in optim(initial_values, fn = mvn_log_likelihood, mvn_log_likelihood_gr,  :
    # non-finite value supplied by optim
    # The fix was to adjust pgtol from the default of 0 to some small number like 1E-12
    # in optim control options

    ll_vals[i] = tryCatch(expr =
                            {
                              -1*stats::optim(initial_values, fn = mvn_log_likelihood_H0,
                                       gr = mvn_log_likelihood_gr_H0,
                                       sample = mvn_sample,
                                       inv_cov_mat = inv_cov_mat,
                                       affected_pollutants = affected_pollutants,
                                       method = "L-BFGS-B",
                                       lower = optimization_limits_list[[i]]$lower,
                                       upper = optimization_limits_list[[i]]$upper)$value
                            },
                          error = function(err)
                          {
                            -1*stats::optim(initial_values, fn = mvn_log_likelihood_H0,
                                     gr = mvn_log_likelihood_gr_H0,
                                     sample = mvn_sample,
                                     inv_cov_mat = inv_cov_mat,
                                     affected_pollutants = affected_pollutants,
                                     method = "L-BFGS-B",
                                     lower = optimization_limits_list[[i]]$lower,
                                     upper = optimization_limits_list[[i]]$upper,
                                     control = list(pgtol = 1E-12))$value
                          })

    # ll_vals[i] = -1*optim(initial_values, fn = mvn_log_likelihood_H0,
    #                       gr = mvn_log_likelihood_gr_H0,
    #                       sample = mvn_sample,
    #                       inv_cov_mat = inv_cov_mat,
    #                       affected_pollutants = affected_pollutants,
    #                       method = "L-BFGS-B",
    #                       lower = optimization_limits_list[[i]]$lower,
    #                       upper = optimization_limits_list[[i]]$upper)$value


  }

  ll = max(ll_vals)

  return(ll)
}

lrt_optimization_H1 = function(mvn_sample, inv_cov_mat,optimization_limits_list)
{
  ll_vals = numeric(length(optimization_limits_list)) # stores the log likelihood

  for(i in 1: length(optimization_limits_list))
  {

    initial_values = generate_initial_vals(optimization_lower = optimization_limits_list[[i]]$lower,
                                           optimization_upper = optimization_limits_list[[i]]$upper)

    # this try catch is because of the following error:
    # Error in optim(initial_values, fn = mvn_log_likelihood, mvn_log_likelihood_gr,  :
    # non-finite value supplied by optim
    # The fix was to adjust pgtol from the default of 0 to some small number like 1E-12
    # in optim control options

    ll_vals[i] = tryCatch(expr =
    {
      -1*stats::optim(initial_values, fn = mvn_log_likelihood,
               mvn_log_likelihood_gr,
               sample = mvn_sample,
               inv_cov_mat = inv_cov_mat,
               method = "L-BFGS-B",
               lower = optimization_limits_list[[i]]$lower,
               upper = optimization_limits_list[[i]]$upper)$value
    },
    error = function(err)
    {
      -1*stats::optim(initial_values, fn = mvn_log_likelihood,
               mvn_log_likelihood_gr,
               sample = mvn_sample,
               inv_cov_mat = inv_cov_mat,
               method = "L-BFGS-B",
               lower = optimization_limits_list[[i]]$lower,
               upper = optimization_limits_list[[i]]$upper,
               control = list(pgtol = 1E-12))$value
    })

    # ll_vals[i] = -1*optim(initial_values, fn = mvn_log_likelihood,
    #                       mvn_log_likelihood_gr,
    #                       sample = mvn_sample,
    #                       inv_cov_mat = inv_cov_mat,
    #                       method = "L-BFGS-B",
    #                       lower = optimization_limits_list[[i]]$lower,
    #                       upper = optimization_limits_list[[i]]$upper)$value


  }

  ll = max(ll_vals)

  return(ll)
}

mvn_log_likelihood_H0 = function(mean, sample, inv_cov_mat, affected_pollutants) # call likelihood with 0s for group of interest
{
  # insert 0s for mean under H0
  mvn_log_likelihood(mean = mean %>%
                       append(rep(0,length(affected_pollutants)),
                              after = affected_pollutants[1]-1),
                     sample = sample, inv_cov_mat = inv_cov_mat)
}


mvn_log_likelihood = function(mean, sample, inv_cov_mat)
{
  return(crossprod(matrix(sample - mean, ncol = 1), inv_cov_mat) %>%
           tcrossprod(matrix(sample - mean, ncol = 1) %>% t))
}

mvn_log_likelihood_gr_H0 = function(B,sample, inv_cov_mat, affected_pollutants)
{
  return(mvn_log_likelihood_gr(B %>% append(rep(0,length(affected_pollutants)),
                                            after = affected_pollutants[1]-1), sample, inv_cov_mat)[-affected_pollutants])
}

mvn_log_likelihood_gr = function(B,sample, inv_cov_mat)
{
  # -1 is for optim function minimzation
  return( -1* tcrossprod(inv_cov_mat, matrix(sample - B, nrow = 1)))
}
################################################################################

sig_code = function(pvalues)
{
  sig_code_vec = character(length = length(pvalues))

  for(i in 1: length(pvalues))
  {
    if(pvalues[i] <= .001)
    {
      sig_code_vec[i] = "***"
    }else if(pvalues[i] <= .01)
    {
      sig_code_vec[i] = "**"
    }else if(pvalues[i] <= .05)
    {
      sig_code_vec[i] = "*"
    }else if(pvalues[i] <= .1)
    {
      sig_code_vec[i] = "."
    }
  }
  return(sig_code_vec)
}

format_scientific = function(pvalues, cutoff = 1E-4)
{
  new_pvalues = pvalues
  digits = -log(cutoff, base = 10)
  for(i in 1: length(pvalues))
  {
    if((pvalues[i] <= cutoff) && (pvalues[i] !=0))
    {
      new_pvalues[i] = format(pvalues[i],digits =digits, scientific = T)
    }else if(pvalues[i] == 0)
    {
      new_pvalues[i] = paste("<", cutoff, sep = "")
    }
    else
    {
      new_pvalues[i] = round(pvalues[i],digits)
    }
  }
  return(new_pvalues)
}

# main caller function to perform inference for FGWQSR
perform_inference = function(ll_models, params_logistic_form, vars,cov_mat,
                             zero_threshold_cutoff, num_sims,reparam_reg, cores,
                             verbose)
{
  num_group_effects = vars$mixture %>% length

  num_in_each_group = sapply(vars$mixture, length)

  num_pollutants = num_in_each_group %>% sum

  pollutant_MLEs = params_logistic_form[2:(2+num_pollutants-1)]

  MLE_effects = listify(pollutant_MLEs, num_in_each_group)

  pollutant_cov_mat = cov_mat[2:(2+num_pollutants-1), 2:(2+num_pollutants-1)]

  # the order of ll_models will be full model, num_pollutants SPLRTs,
  # then num_group_effects 1 group LRTs

  # create matrix with group_num and pollutant_num to call in future_apply

  inference_order = matrix(nrow = num_pollutants+num_group_effects,ncol = 2)
  colnames(inference_order) = c("Group Num", "Pollutant Num")

  # populate group number column
  inference_order[1:num_pollutants,1] = rep(1:num_group_effects, times = num_in_each_group)
  inference_order[(num_pollutants+1): nrow(inference_order),1] = 1:num_group_effects

  # populate pollutant numbers
  inference_order[1:num_pollutants,2] = sapply(num_in_each_group, FUN = function(x) 1:x ) %>% unlist
  # group effects are NA

  # get LRT statistics from FGWQSR
  # first ll is full model, all other ll's are sub models
  lrts_FG = -2*(ll_models[-1] -ll_models[1])

  # make sure numerically negative LRTs are set to 0
  lrts_FG = ifelse(lrts_FG < 0,0, lrts_FG)

  # perform inference on single pollutant and group effects
  lrt_dists = vector(mode = "list", length = nrow(inference_order))

  # set a threshold for when LRT simulated distributions will not be computed
  LRT_threshold = 1E-8

  # run models in parallel, verbose option

  if(verbose == TRUE)
  {

    p <- progressr::progressor(steps = nrow(inference_order))
    progressr::handlers("progress")
  }

  for(i in 1: nrow(inference_order))
  {
    pollutant_num = inference_order[i,2]
    null_effects = MLE_effects

    # if a LRT is less than a threshold, do not waste time computing LRT
    # distribution
    if(lrts_FG[i] <  LRT_threshold)
    {
      lrt_dists[[i]] = Inf

    }else if(!is.na(pollutant_num)) # if SPLRT
    {
      group_num = inference_order[i,1]
      null_effects[[group_num]][pollutant_num] = 0
      lrt_dists[[i]] = mvn_lrt_simulation(effects = null_effects, cov_mat = pollutant_cov_mat,
                                          group_num = inference_order[i,1],
                                          pollutant_num =  pollutant_num,
                                          num_sims = num_sims,
                                          zero_threshold_cutoff = zero_threshold_cutoff,
                                          cores = cores)

    }else  # if one group LRT, simulate with group effects being null
    {
      group_num = inference_order[i,1]
      null_effects[[group_num]] = rep(0, num_in_each_group[group_num])
      lrt_dists[[i]] = mvn_lrt_simulation(effects = null_effects, cov_mat = pollutant_cov_mat,
                                          group_num = inference_order[i,1],
                                          pollutant_num =  pollutant_num,
                                          num_sims = num_sims,
                                          zero_threshold_cutoff = zero_threshold_cutoff,
                                          cores = cores)
    }

    if(verbose == TRUE)
    {
      p()
    }
  }

  # make sure values in lrt_dists are not numerically negative, make them 0
  lrt_dists = lapply(lrt_dists, FUN = function(x) ifelse(x != Inf & (x < 0), 0, x))

  # get pvalues for SPLRTs and 1 group LRTs

  pvals_polls = sapply(1:nrow(inference_order), FUN = function(i)
    {
    pvalue = 0
    if(lrt_dists[[i]][1] == Inf)
    {
      pvalue = 1
    }else
    {
      pvalue = sum( lrts_FG[i] <= lrt_dists[[i]]) / num_sims
    }
    return( pvalue )
  })

  # perform inference on intercept and confounders

  # subset the parameters and the variances from observed FI matrix
  adj_params = params_logistic_form[-c(2:(2+num_pollutants-1))]

  adj_params_var = cov_mat[-c(2:(2+num_pollutants-1)),
                               -c(2:(2+num_pollutants-1))] %>% as.matrix %>% diag

  # compute Z
  Z_adj_params = adj_params / sqrt(adj_params_var)

  p_val_adj = (2*(1-stats::pnorm(abs(Z_adj_params))))

  # compute CIs

  ci_95_lower = adj_params - stats::qnorm(.975)*sqrt(adj_params_var)
  ci_95_upper = adj_params + stats::qnorm(.975)*sqrt(adj_params_var)

  group_index_frame = data.frame(Estimate = reparam_reg$group_index,
                           LRT = lrts_FG[(num_pollutants+1): nrow(inference_order)],
                           "P-value" = pvals_polls[(num_pollutants+1): nrow(inference_order)],
                           " " = sig_code(pvals_polls[(num_pollutants+1):
                                                        nrow(inference_order)]%>% as.numeric),
                           check.names = F)

  rownames(group_index_frame) = paste("Mixture Effect", 1:num_group_effects)

  # weights / single pollutant LRTs

  weight_frame = data.frame("Weight Estimate" = reparam_reg$weights %>% unlist,
                            LRT = lrts_FG[1:num_pollutants],
                            "P-value" = pvals_polls[1:num_pollutants],
                            " " = sig_code(pvals_polls[1:num_pollutants]%>% as.numeric),
                            check.names = F)

  rownames(weight_frame) = names(pollutant_MLEs)[1:num_pollutants]

  adj_param_frame = data.frame(Estimate = adj_params,
                               SE = sqrt(adj_params_var),
                               Z = Z_adj_params,
                               "P(Z > |z|)" = p_val_adj,
                               ci_95_lower = ci_95_lower,
                               ci_95_upper = ci_95_upper,
                               " " = sig_code(p_val_adj %>% as.numeric),
                               check.names = F)


  return(list(group_index_frame = group_index_frame,
              weight_frame = weight_frame,
              adj_param_frame = adj_param_frame))
}









