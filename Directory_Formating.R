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

## Now we need to add some "dummy" text files in each subdir so we can add these to github for
## the team to see the structure.

mod.paths <- list.dirs("./Modules/")
target_paths <- mod.paths[grep("(Code|Data|Presentations)$", mod.paths)]

# Create .gitkeep files in those directories
for(path in target_paths) {
  file.create(file.path(path, "dummy.txt"))
  sometext <- "Dummy file. Delete when ready."
  writeLines(sometext, file.path(path, "dummy.txt"))
}
