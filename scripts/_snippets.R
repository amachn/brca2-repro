# time-to-run timing snippet:
ST <- Sys.time()
# insert code to be timed here
ET <- Sys.time()
time_taken <- paste("expr_list reduction:", difftime(ET, ST, units = "secs"), "seconds")

# base R environment reset snippet:
ip <- as.data.frame(installed.packages())
ip <- ip[!(ip[, "Priority"] %in% c("base", "recommended")), ]
path.lib <- unique(ip$LibPath)
pkgs.to.remove <- ip[, 1]
sapply(pkgs.to.remove, remove.packages, lib = path.lib)
