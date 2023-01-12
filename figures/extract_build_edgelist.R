r <- readRDS("../resources/remembr.RDS")

write.network.files <- function (case.name) {
    edgelist <- r[[paste0(case.name, ".ortho.edgelist.90")]][["long.run"]]
    colnames(edgelist) <- c("target", "source", "ignore", "eco.type")

    ## remove the loops
    edgelist <- edgelist[edgelist$source != edgelist$target,]

    nodelist <- unique(r[[paste0(case.name, ".fit.df")]]$subreddit)
    nodelist <- data.frame(id=nodelist,
                            label=nodelist)
    ## nodelist <- data.frame(id=unique(c(edgelist$source, edgelist$target)),
    ##                        label=unique(c(edgelist$source, edgelist$target)))

    write.table(edgelist, file=paste0(case.name, "_edgelist.csv"), row.names=FALSE, sep="\t")
    write.table(nodelist, file=paste0(case.name, "_nodelist.csv"), row.names=FALSE, sep="\t")
}

write.network.files("seattle")
write.network.files("wallpaper")
write.network.files("design")


