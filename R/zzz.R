# Welcome information of the package
.onAttach <- function(libname, pkgname) {
  version = utils::packageDescription(pkgname, fields = "Version")
  msg = paste0(
    "========================================
",
    pkgname,
    " version ",
    version,
    "
Github page: https://github.com/Ronlee12355/HNSCclassifier
Report bugs: https://github.com/Ronlee12355/HNSCclassifier/issues

This message can be suppressed by:
  suppressPackageStartupMessages(library(HNSCclassifier))
========================================
"
  )

  packageStartupMessage(msg)
}
