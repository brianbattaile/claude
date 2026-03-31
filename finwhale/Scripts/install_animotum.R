repos <- c("https://ianjonsen.r-universe.dev", "https://cloud.r-project.org")

pkgs <- c("geographiclib", "traipse", "aniMotum")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = repos, dependencies = TRUE)
  } else {
    cat(pkg, "already installed (", as.character(packageVersion(pkg)), ")\n")
  }
}

library(aniMotum)
cat("aniMotum version:", as.character(packageVersion("aniMotum")), "\n")
