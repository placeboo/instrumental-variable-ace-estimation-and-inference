## input is hmat
## M = (l,b,-A)
## cols: ineq(l=0)/eq(l=1); right-hand-side;p(x0,y0|z1) p(x0,y1|z1) p(x1,y0|z1) ... p(x1,y1 | zp)
## Ax<=b; Ax=b
## rows: each surface

descrip_hmat = function(mat){
      N_surface = nrow(mat)
      N_p = ncol(mat)

      eq_index = which(mat[,1] == 1)
      ineq_index = which(mat[,1] == 0)

      ineq_mat = mat[ineq_index, ]
      eq_mat = mat[eq_index, ]

      ## inequalities
      ineq_1 = ineq_mat[1, ]
      ## each type is organized as a vector: b, #-1, #0, #1
      tmp_vec = c("b"=0, "-1"=0, "0" = 0, "1"= 0)
      type_1 = c("b" = ineq_1[2],table(ineq_1[3: N_p]))
      tmp_vec[names(type_1)] = type_1
      type_1 = tmp_vec

      type_list = list(type_1)
      type_num = c(1) # record the number of each type
      type_index_list = list(1) # record the surface index for each type

      type_index = 1 # number of different types
      for(i in 2: nrow(ineq_mat)){
            tmp_vec = c("b"=0, "-1"=0, "0" = 0, "1"= 0) # c(b, #-1, #0, #1)
            ineq_i = ineq_mat[i, ]
            type_i = c("b"=ineq_i[2],table(ineq_i[3: N_p]))# c(b, #-1, #0, #1)
            tmp_vec[names(type_i)] = type_i

            type_i = tmp_vec
            # check if the type_i in the type list
            in_list = type_list %in% list(type_i) # a vecotr of T/F
            if(any(in_list)){ # type_i in the list
                  # add 1 to the certain old type
                  type_num[which(in_list)] = 1 + type_num[which(in_list)]
                  type_index_list[[which(in_list)]] = c(type_index_list[[which(in_list)]], i)
            }else{ # new type appear
                  type_index = type_index + 1
                  # add the type into type list
                  type_list[[type_index]] = type_i
                  # add 1 to the new type number
                  type_num[type_index] = 1
                  type_index_list[[type_index]] = i
            }
      }

      type_mat = do.call(rbind, type_list)
      type_mat = cbind(type_mat, num = type_num)

      result_list = list(eq_num = length(eq_index), ineq = type_mat, ineq_type_index = type_index_list)
      return(result_list)
}



