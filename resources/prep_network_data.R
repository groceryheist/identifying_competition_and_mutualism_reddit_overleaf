print(getwd())
source("resources/preamble.R")
prep.network <- function(irf, name){

    nodes <- unique(irf,by=c("impulse.subreddit"))[,.(impulse.subreddit)]
    nodes[,":="(impulse.subreddit = gsub("^X","",impulse.subreddit))]
    nodes[,":="(id = impulse.subreddit, label=impulse.subreddit)]
    print(nodes)

    fwrite(nodes,paste0(name,'_node_list','.csv'))

    irf <- plot.irf.ols.data(irf)
    
    irf <- irf[irf.lower < 0,'eco.type':='competition']
    irf <- irf[irf.upper > 0,'eco.type':='mutualism']
    irf <- irf[!(is.na(eco.type)),.(eco.type=first(eco.type)),by=c("impulse.subreddit","response.subreddit")]

    edge.list <- unique(irf[!is.na(eco.type)] ,by=c("impulse.subreddit", "response.subreddit","eco.type"))[,.(impulse.subreddit,response.subreddit,'ignore'="",eco.type)]

    edge.list <- edge.list[,":="(impulse.subreddit = gsub("^r.","",impulse.subreddit),
                    response.subreddit = gsub("^r.","",response.subreddit))]
    setnames(edge.list, old=c("impulse.subreddit","response.subreddit"),new=c('source','target'))


    fwrite(edge.list,paste0(name,'_edge_list','.csv'))
}

prep.network(law.var$irf.data.author_cluster_101_tf, "law")
prep.network(redpill.var$irf.data.author_cluster_320_tf, "redpill")
prep.network(mental.var$irf.data.author_cluster_468_tf, "mental")
prep.network(cod.var$irf.data.author_cluster_599_tf, "cod")
prep.network(data.var$irf.data.author_cluster_631_tf, "data")
prep.network(realestate.var$irf.data.author_cluster_575_tf, "realestate")
prep.network(watches.var$irf.data.author_cluster_513_tf, "watches")

## prep.network(gun.var$irf.data.author_cluster_829_tf, "gun")
## prep.network(diet.var$irf.data.author_cluster_1_tf, "diet")
## prep.network(art.var$irf.data.author_cluster_1636_tf, "art")
## prep.network(toys.var$irf.data.author_cluster_23_tf, "toy")
## prep.network(energy.var$irf.data.author_cluster_1524_tf, "energy")
## prep.network(tanks.var$irf.data.author_cluster_1763_tf, "tank")
## prep.network(tumblr.var$irf.data.author_cluster_1878_tf, "tumblr")
