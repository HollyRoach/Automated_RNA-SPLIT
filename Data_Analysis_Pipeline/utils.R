

#create function to replace choose.dir so that it is compatible with mac OS
choose_dir <- function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}