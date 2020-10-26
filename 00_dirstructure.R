
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
inspectdir <- file.path(resdir, "data_inspection")
moddir <- file.path(resdir, "data_modeling")

if (dir.exists(inspectdir)) {
  warning(paste(inspectdir, 'already exists'))
} else {
  warning(paste('Create new directory:',inspectdir))
  dir.create(inspectdir)
}


if (dir.exists(moddir)) {
  warning(paste(moddir, 'already exists'))
} else {
  warning(paste('Create new directory:',moddir))
  dir.create(moddir)
}