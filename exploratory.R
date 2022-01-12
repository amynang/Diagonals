temp1 = tempfile()
download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Soil+food+webs%2C+Successions+-+Jacobian+Matrices?entryid=1acea44e-44da-4225-8e98-203d742a0a82&output=zip.zip",temp1)

temp2 = tempfile()
download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Antarctic+food+webs+-+Jacobian+matrices?entryid=1d1e4253-552a-421e-9bdc-20d70f683cc5&output=zip.zip",temp2)

# temp3 = tempfile()
# download.file("https://ramadda.data.bas.ac.uk/repository/entry/show/Polar+Data+Centre/DOI/Soil+food+webs%2C+agro+and+native+-+Jacobian+matrices?entryid=7b6dd454-d248-4264-b8a0-2de07a8dc882&output=zip.zip",temp3)



jacob = vector(mode = "list", 35)
for (i in 1:33) { 
jacob[[i]] = read.csv(unzip(temp1)[i])
rownames(jacob[[i]]) = jacob[[i]][,1]
jacob[[i]] = as.matrix(jacob[[i]][,-1])
colnames(jacob[[i]]) = rownames(jacob[[i]])
}
#rename rows &remove first column!!!!!!!!
jacob[[34]] = read.csv(unzip(temp2)[1])
jacob[[35]] = read.csv(unzip(temp2)[3])

for (i in 34:35) {  
rownames(jacob[[i]]) = jacob[[i]][,1]
jacob[[i]] = as.matrix(jacob[[i]][,-1])
colnames(jacob[[i]]) = rownames(jacob[[i]])
}

jacob = jacob[-17]

saveRDS(jacob, file="jacob.RData")
readRDS("jacob.RData")


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
b = max(unstable.jacob[[30]][unstable.jacob[[30]]<0])

unstable.jacob[[30]][1,1] = -0.005660854
max(Re(eigen(unstable.jacob[[30]])$values))
diag(unstable.jacob[[30]]) = 0
max(Re(eigen(unstable.jacob[[30]])$values))

View(unstable.jacob[[30]])

m = nrow(unstable.jacob[[30]])

selfreg = data.frame(matrix(ncol = m+1, nrow = 100))

max(unstable.jacob[[30]][unstable.jacob[[30]]<0])
selfreg[,1] = seq(from = 0, to = -0.005660854*9900, by = -0.005660854*100)

for (i in 1:100) { 
  
  diag(unstable.jacob[[30]]) = 0
  
  for (j in 1:m) { 
    diag(unstable.jacob[[30]]) = 0
    unstable.jacob[[30]][j,j] = selfreg[i,1]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
    selfreg[i,1+j] = max(Re(eigen(unstable.jacob[[30]])$values)) 
    diag(unstable.jacob[[30]]) = 0
  }
}

plot(-selfreg[,1],selfreg[,2],type="l")
plot(-selfreg[,1],selfreg[,3],type="l")
plot(-selfreg[,1],selfreg[,4],type="l")
plot(-selfreg[,1],selfreg[,5],type="l")
plot(-selfreg[,1],selfreg[,6],type="l")
plot(-selfreg[,1],selfreg[,7],type="l")
plot(-selfreg[,1],selfreg[,8],type="l")
plot(-selfreg[,1],selfreg[,9],type="l")






m = nrow(unstable.jacob[[34]])

selfreg = data.frame(matrix(ncol = m+1, nrow = 100))

max(unstable.jacob[[34]][unstable.jacob[[34]]<0])
selfreg[,1] = seq(from = 0, 
                  to = max(unstable.jacob[[34]][unstable.jacob[[34]]<0])*99, 
                  by = max(unstable.jacob[[34]][unstable.jacob[[34]]<0])*1)

for (i in 1:100) { 
  
  diag(unstable.jacob[[34]]) = 0
  
  for (j in 1:m) { 
    diag(unstable.jacob[[34]]) = 0
    unstable.jacob[[34]][j,j] = selfreg[i,1]#runif(1, min=a, max = b)#diagonals[[72]][9,i]
    selfreg[i,1+j] = max(Re(eigen(unstable.jacob[[34]])$values)) 
    diag(unstable.jacob[[34]]) = 0
  }
}

plot(-selfreg[,1],selfreg[,2])
plot(-selfreg[,1],selfreg[,3])
plot(-selfreg[,1],selfreg[,4])
plot(-selfreg[,1],selfreg[,5])
plot(-selfreg[,1],selfreg[,6])
plot(-selfreg[,1],selfreg[,7])
plot(-selfreg[,1],selfreg[,8])
plot(-selfreg[,1],selfreg[,9])

plot(-selfreg[,1],selfreg[,12])
plot(-selfreg[,1],selfreg[,13])
plot(-selfreg[,1],selfreg[,14])
plot(-selfreg[,1],selfreg[,15])
plot(-selfreg[,1],selfreg[,16])
plot(-selfreg[,1],selfreg[,17])
plot(-selfreg[,1],selfreg[,18])
plot(-selfreg[,1],selfreg[,19])








library(beepr)

random.move = function(diagonal.element, solution) {
  #select a diagonal element to change
  sel = sample(diagonal.element, 1)
  #keep the *original maximum (min) self regulation
  min = jacob[[30]][sel,sel]
  #min = -6
  #draw a random value from the [min,0] range
  new = runif(1, min = min, max = 0)
  #replace that element with the random value
  solution[sel,sel] = new
  new.solution = solution
  return(new.solution)
}


compute.energy = function(jacobian){
  
  stab = max(Re(eigen(jacobian)$values))
  return(stab)
}

simmulated_annealing = function(jacobian, t.init = 1, t.decrease = 0.999, t.final = .999^100000){
  
  
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
    new.sol = random.move(diagonal.element, solution)
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
    #diagonals[nrow(diagonals) + 1,] = diag(solution)
    #diagonals = rbind(diagonals, diag(solution))
    #print(diagonals)
  }# end of the loop on t
  
  return(solution.best)
}
stabler = simmulated_annealing(unstable.jacob[[30]]) ; beep(9)


stabler = simmulated_annealing(unstable.jacob[[30]]) ; beep(9)
