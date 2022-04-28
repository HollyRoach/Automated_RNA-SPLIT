#' Replace choose.dir so that it is compatible with mac OS
#' 
#' @param caption A text string to display as the folder caption
#' 
#' @return A string of the file path
choose_dir <- function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

#' Remove all objects from the environment except a subset passed in
#' 
#' @param except A list of objects as characters
clean_env <- function(except) {
  parent = sys.frame(-1)
  all_objs <- ls(envir = parent)
  objs_to_remove <- all_objs[!(all_objs %in% except)]
  rm(list = objs_to_remove, envir = parent)
}