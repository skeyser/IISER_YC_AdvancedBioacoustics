## Formatting the Rproject
subdirs <- c("Code", "Data", "Presentations")

paths <- list.dirs("./Modules/", recursive = F)

for(i in 1:length(paths)){
  path.tmp <- paths[i]
  if(all(dir.exists(paste0(path.tmp, "/", subdirs)))){
    next
  } else {
    if(!dir.exists(paste0(path.tmp, "/", subdirs[1]))){
      dir.create(paste0(path.tmp, "/", subdirs[1]))
    }
    if(!dir.exists(paste0(path.tmp, "/", subdirs[2]))){
      dir.create(paste0(path.tmp, "/", subdirs[2]))
    }
    if(!dir.exists(paste0(path.tmp, "/", subdirs[3]))){
      dir.create(paste0(path.tmp, "/", subdirs[3]))
    }
  }
}
