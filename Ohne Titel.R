source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

#install.packages("quadprog")
#install.packages("scdensity")
#install.packages("LaplacesDemon")
#install.packages("devtools")
#devtools::install_github("ccombesGG4/AmoRosoDistrib")


dat <- rnorm(100)

densityMclust(dat, plot = F)

?densityMclust


safe_execute <- function(expr, object_name, data_vector) {
  tryCatch(
    {
      result <- eval(bquote(.(expr)), envir = list(dat = data_vector))
      return(result)
    },
    error = function(e) {
      cat(paste("Error with fitting", object_name, ":", e$message, ";\n",
                "Other methods were still fit.\n"))
      return(NA) # assigns NA if method failed
    }
  )
}
mnorm <- safe_execute(quote(
  densityMclust(dat,plot=F)), "mnorm", dat)

              