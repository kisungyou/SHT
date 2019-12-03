.pkgenv <- new.env(parent = emptyenv())

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = packageVersion("SHT")
  
  ## Print on Screen
  packageStartupMessage("**------------------------------------------------------**")
  packageStartupMessage("** MMMMMMMMMMMMMMMMMMMMM MMMMMMMMMMMM MMMMMMMMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMMP       MM M  MMMMM  MM M        MMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMM  mmmmm..M M  MMMMM  MM MMMM  MMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMM.....   YM M         MM MMMM  MMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMMMMMMMM.  M M  MMMMM  MM MMMM  MMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMM. .MMM   M M  MMMMM  MM MMMM  MMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMMb.     .dM M  MMMMM  MM MMMM  MMMMMMMMMMMMM ")
  packageStartupMessage("** MMMMMMMMMMMMMMMMMMMMM MMMMMMMMMMMM MMMMMMMMMMMMMMMMMMM ")
  packageStartupMessage("** ")
  packageStartupMessage("**         Statistical Hypothesis Testing Toolbox")
  packageStartupMessage("**")
  packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Kisung You (kyoustat@gmail.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("**------------------------------------------------------**")
}

.onUnload <- function(libpath) {
  library.dynam.unload("SHT", libpath)
}