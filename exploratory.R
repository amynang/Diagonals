library(ggplot2)
library(tidyverse)

# temp1 = tempfile()
# download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Soil+food+webs%2C+Successions+-+Jacobian+Matrices?entryid=1acea44e-44da-4225-8e98-203d742a0a82&output=zip.zip",temp1)
# 
# temp2 = tempfile()
# download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Antarctic+food+webs+-+Jacobian+matrices?entryid=1d1e4253-552a-421e-9bdc-20d70f683cc5&output=zip.zip",temp2)
# 
# temp3 = tempfile()
# download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Soil+food+webs%2C+agro+and+native+-+Jacobian+matrices?entryid=7b6dd454-d248-4264-b8a0-2de07a8dc882&output=zip.zip",temp3)
# 
# colnames(jacob[[1]]) = case_when()
# 
# jacob = vector(mode = "list", 35)
# for (i in 1:33) { 
# jacob[[i]] = read.csv(unzip(temp1)[i])
# jacob[[i]][,1] = case_when(jacob[[i]][,1] == "PRMI" ~ "Predatory mites",
#                            jacob[[i]][,1] == "PRCO" ~ "Predatory collembolans",
#                            jacob[[i]][,1] == "NEMI" ~ "Nematode feeding mites",
#                            jacob[[i]][,1] == "PRNM" ~ "Predatory nematodes",
#                            jacob[[i]][,1] == "AMOE" ~ "Amoebae",
#                            jacob[[i]][,1] == "COLL" ~ "Collembolans",
#                            jacob[[i]][,1] == "CRYP" ~ "Cryptostigmatic mites",
#                            jacob[[i]][,1] == "NCRY" ~ "Noncryptostigmatic mites",
#                            jacob[[i]][,1] == "FUNE" ~ "Fungivorous nematodes",
#                            jacob[[i]][,1] == "FLAG" ~ "Flagellates",
#                            jacob[[i]][,1] == "BANE" ~ "Bacteriophagous nematodes",
#                            jacob[[i]][,1] == "BAMI" ~ "Bacteriophagous mites",
#                            jacob[[i]][,1] == "PHNE" ~ "Phytophageous nematodes",
#                            jacob[[i]][,1] == "FUNG" ~ "Saprophytic fungi",
#                            jacob[[i]][,1] == "BACT" ~ "Bacteria",
#                            jacob[[i]][,1] == "ROOT" ~ "Roots",
#                            jacob[[i]][,1] == "DETR" ~ "Detritus")
# rownames(jacob[[i]]) = jacob[[i]][,1]
# jacob[[i]] = as.matrix(jacob[[i]][,-1])
# colnames(jacob[[i]]) = rownames(jacob[[i]])
# }
# #rename rows &remove first column!!!!!!!!
# jacob[[34]] = read.csv(unzip(temp2)[1])
# jacob[[35]] = read.csv(unzip(temp2)[3])
# 
# for (i in 34:35) {  
# rownames(jacob[[i]]) = jacob[[i]][,1]
# jacob[[i]] = as.matrix(jacob[[i]][,-1])
# colnames(jacob[[i]]) = rownames(jacob[[i]])
# }
# 
# jacob = jacob[-17]
# 
# saveRDS(jacob, file="jacob.RData")
jacob = readRDS("jacob.RData")


for (i in 1:34) { 
#diag(jacob[[i]]) = 0
lambda = max(Re(eigen(jacob[[i]])$values))
  print(lambda)
}

unstable.jacob = jacob

for (i in 1:34) { 
  diag(unstable.jacob[[i]]) = 0
  lambda = max(Re(eigen(unstable.jacob[[i]])$values))
  print(lambda)
}
# b = max(unstable.jacob[[30]][unstable.jacob[[30]]<0])
# 
# unstable.jacob[[30]][1,1] = -0.005660854
# max(Re(eigen(unstable.jacob[[30]])$values))
# diag(unstable.jacob[[30]]) = 0
# max(Re(eigen(unstable.jacob[[30]])$values))
# 
# View(unstable.jacob[[30]])
# 
# m = nrow(unstable.jacob[[30]])
# 
# selfreg = data.frame(matrix(ncol = m+1, nrow = 100))
# 
# max(unstable.jacob[[30]][unstable.jacob[[30]]<0])
# 
# selfreg[,1] = seq(from = 0, to = -0.005660854*9900, by = -0.005660854*100)
# 
# for (i in 1:100) { 
#   
#   diag(unstable.jacob[[30]]) = 0
#   
#   for (j in 1:m) { 
#     diag(unstable.jacob[[30]]) = 0
#     unstable.jacob[[30]][j,j] = selfreg[i,1]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
#     selfreg[i,1+j] = max(Re(eigen(unstable.jacob[[30]])$values)) 
#     diag(unstable.jacob[[30]]) = 0
#   }
# }
# 
# plot(-selfreg[,1],selfreg[,2],type="b")
# plot(-selfreg[,1],selfreg[,3],type="b")
# plot(-selfreg[,1],selfreg[,4],type="b")
# plot(-selfreg[,1],selfreg[,5],type="b")
# plot(-selfreg[,1],selfreg[,6],type="b")
# plot(-selfreg[,1],selfreg[,7],type="b")
# plot(-selfreg[,1],selfreg[,8],type="b")
# plot(-selfreg[,1],selfreg[,9],type="b")

# starting from zero diagonals going to max values for individual nodes
for (m in 1:34) {  
  png(paste0("unstable", m ,".png"), 800, 900, pointsize = 18)
  par(mfrow=c(ceiling(sqrt(nrow(unstable.jacob[[m]]))),
              ceiling(sqrt(nrow(unstable.jacob[[m]])))))
  
  selfreg = vector(mode="list", length = nrow(unstable.jacob[[m]]))
  names(selfreg) = colnames(unstable.jacob[[m]])
  for (i in 1:nrow(unstable.jacob[[m]])) { 
    
    selfreg[[i]] = data.frame(matrix(ncol = 2, nrow = 100))
    colnames(selfreg[[i]]) = c("SelfRegulation","lambda")
    diag(unstable.jacob[[m]]) = 0
    
    
    for (j in 1:100) { 
      diag(unstable.jacob[[m]]) = 0
      selfreg[[i]][ ,"SelfRegulation"] = seq(from = 0, to = jacob[[m]][i,i], length.out = 100)
      unstable.jacob[[m]][i,i] = selfreg[[i]][j,"SelfRegulation"]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
      selfreg[[i]][j,"lambda"] = max(Re(eigen(unstable.jacob[[m]])$values)) 
      diag(unstable.jacob[[m]]) = 0
    }
    
    plot(-selfreg[[i]]$SelfRegulation,
         selfreg[[i]]$lambda,
         type="b",
         main = names(selfreg)[i]
         ,xlab = "Self Regulation"
         ,ylab = "Lambda"
    )
  }
  title(m, line = -1, side = 3, outer = TRUE)
  
  #png(paste0("unstable", m ,".png"))
  #sjp.frq(data[,i],title = names(data)[i])
  dev.off()
}

# starting from max diagonals going to zero for individual nodes
jacob.2 = jacob
for (m in 1:34) { 
  png(paste0("stable", m ,".png"), 800, 900, pointsize = 18)
  par(mfrow=c(ceiling(sqrt(nrow(jacob.2[[m]]))),
              ceiling(sqrt(nrow(jacob.2[[m]])))))
  
  selfreg = vector(mode="list", length = nrow(jacob.2[[m]]))
  names(selfreg) = colnames(jacob.2[[m]])
  for (i in 1:nrow(jacob.2[[m]])) { 
    
    selfreg[[i]] = data.frame(matrix(ncol = 2, nrow = 100))
    colnames(selfreg[[i]]) = c("SelfRegulation","lambda")
    diag(jacob.2[[m]]) = diag(jacob[[m]])
    
    
    for (j in 1:100) { 
      diag(jacob.2[[m]]) = diag(jacob[[m]])
      selfreg[[i]][ ,"SelfRegulation"] = seq(from = jacob[[m]][i,i], to = 0, length.out = 100)
      jacob.2[[m]][i,i] = selfreg[[i]][j,"SelfRegulation"]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
      selfreg[[i]][j,"lambda"] = max(Re(eigen(jacob.2[[m]])$values)) 
      #diag(jacob.2[[m]]) = diag(jacob[[m]])
    }
    
    plot(-selfreg[[i]]$SelfRegulation,
         -selfreg[[i]]$lambda,
         type="b",
         main = names(selfreg)[i]
         ,xlab = "Self Regulation"
         ,ylab = "-Lambda"
    )
  }
  title(m, line = -1, side = 3, outer = TRUE)
  dev.off()
}

# p1 = ggplot(selfreg[[1]]) +
#   theme_bw() +
#   geom_point(aes(-SelfRegulation, lambda)) +
#   ggtitle(names(selfreg)[1])
# 
# 
# 
# 
# 
# m = nrow(unstable.jacob[[34]])
# 
# selfreg = data.frame(matrix(ncol = m+1, nrow = 100))
# 
# max(unstable.jacob[[34]][unstable.jacob[[34]]<0])
# selfreg[,1] = seq(from = 0, 
#                   to = max(unstable.jacob[[34]][unstable.jacob[[34]]<0])*99, 
#                   by = max(unstable.jacob[[34]][unstable.jacob[[34]]<0])*1)
# 
# for (i in 1:100) { 
#   
#   diag(unstable.jacob[[34]]) = 0
#   
#   for (j in 1:m) { 
#     diag(unstable.jacob[[34]]) = 0
#     unstable.jacob[[34]][j,j] = selfreg[i,1]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
#     selfreg[i,1+j] = max(Re(eigen(unstable.jacob[[34]])$values)) 
#     diag(unstable.jacob[[34]]) = 0
#   }
# }
# 
# plot(-selfreg[,1],selfreg[,2])
# plot(-selfreg[,1],selfreg[,3])
# plot(-selfreg[,1],selfreg[,4])
# plot(-selfreg[,1],selfreg[,5])
# plot(-selfreg[,1],selfreg[,6])
# plot(-selfreg[,1],selfreg[,7])
# plot(-selfreg[,1],selfreg[,8])
# plot(-selfreg[,1],selfreg[,9])
# 
# plot(-selfreg[,1],selfreg[,12])
# plot(-selfreg[,1],selfreg[,18])
# plot(-selfreg[,1],selfreg[,14])
# plot(-selfreg[,1],selfreg[,15])
# plot(-selfreg[,1],selfreg[,16])
# plot(-selfreg[,1],selfreg[,17])
# plot(-selfreg[,1],selfreg[,18])
# plot(-selfreg[,1],selfreg[,19])
# 




thelot = NULL


library(beepr)

compute.energy = function(jacobian){
  
  stab = max(Re(eigen(jacobian)$values))
  return(stab)
}


# N <- 1e6  # total number of rows to preallocate--possibly an overestimate
# 
# thelot = data.frame(matrix(NA,    # Create empty data frame
#                            nrow = N,
#                            ncol = 8))
# colnames(thelot) = colnames(jacob[[30]])
# thelot$Stability = NA
# thelot[1, ] = c(diag(jacob[[30]]), max(Re(eigen(jacob[[30]])$values)))
k = 1
diagonals = t(data.frame(diag(jacob[[k]])))
rownames(diagonals) = NULL

simmulated_annealing = function(jacobian, t.init = 1, t.decrease = 0.999, t.final = .999^100000){
  #thelot[1, ] = c(diag(jacobian), max(Re(eigen(jacobian)$values)))
  
  #counter = 2 
  
  diagonal.element = 1:length(diag(jacobian))
  #nb_l = 1:length(which(jacobian!=0))
  solution = jacobian
  energy = compute.energy(solution)
  # these two things are used to keep in memory the best solutions found
  energy.best = energy
  solution.best = solution
  t = t.init
  while(t>t.final){
    # 1) decrease temperature
    t = t*t.decrease
    
    # 2) make a swap between two elements
    
    #select a diagonal element to change
    sel = sample(diagonal.element, 1)
    #keep the *original maximum (min) self regulation
    min = jacobian[sel, sel]
    #min = -6
    #draw a random value from the [min,0] range
    new = runif(1, min = min, max = 0)
    #replace that element with the random value
    solution[sel,sel] = new
    new.sol = solution
    #new.sol = random.move(diagonal.element, solution)
    new.energy = compute.energy(new.sol)
    
    # 3) accept or reject the move with metropolis hasting criteria
    if (new.energy < energy){
      # here you made a good move, so new.energy and new.solution 
      # are better than the old ones => replace by the new ones
      energy = new.energy
      solution  = new.sol
      if (energy < energy.best){
        # lucky you, you found a new best solution
        # save it for later
        energy.best = energy
        solution.best = solution
      }
    }else{
      # you did not increased the energy
      # still you want to be able to accept the solution
      # but with a given probability that decreases with t
      if (exp((energy - new.energy)/t) > runif(1)){
        # here the move is a bad one, but I still accept it
        energy = new.energy
        solution = new.sol
      }
    }# end of 3 
  
    diagonals = rbind(diagonals, diag(solution.best))
    
  }# end of the loop on t
  print(max(Re(eigen(jacobian)$values)))
  print(max(Re(eigen(solution.best)$values)))
  return(diagonals)
}
stabler = simmulated_annealing(jacob[[k]]) ; beep(9)
#jacob[[30]];stabler

jacob.3 = jacob[[k]]
Stability = vector(mode = "numeric", 100001)
for (i in 1:100001) {
  diag(jacob.3) = stabler[i,]
  Stability[i] = max(Re(eigen(jacob.3)$values))
}
stablerr = cbind(stabler, Stability)
plot(stablerr[,9])
