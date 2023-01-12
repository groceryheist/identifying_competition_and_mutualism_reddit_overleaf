library(lmtest)
#library(dotwhisker)
library(scales)
#library(stargazer)
library(xtable)
library(data.table)
library(knitr)
library(gtable)

format.percent <- function(x,digits=1) {paste(round(x*100,digits),"\\%",sep='')}

## if(sessionInfo()$R.version$major != '4'){
##     suppressMessages(library("lemon",lib.loc="resources"))
## #    suppressMessages(library("ggplot2",lib.loc="resources"))
## } else {
##     suppressMessages(library("lemon"))
##     library(ggplot2)
## }    

library(ggplot2)

if(length(opts_chunk$get("dev")) > 0){
    if(opts_chunk$get("dev") == 'tikz'){
        base_size <- 12
        theme_set(theme_bw(base_size=base_size,base_line_size=base_size/18, base_rect_size=base_size/18))

    } else {
        theme_set(theme_bw(base_size=8))
    }
}
#theme_set(theme_bw(base_size=base_size,base_line_size=base_size/18, base_rect_size=base_size/18))

## r <- readRDS("resources/remember.RDS")
## r2 <- readRDS("resources/remembr.RDS")
## attach(r,warn.conflicts=F)
## attach(r2,warn.conflicts=F)
r <- readRDS("resources/var_ols_remember.RDS")
attach(r,warn.conflicts=F)
rop <- readRDS("resources/remember_overlaps.RDS")
attach(rop,warn.conflicts=F)

law.var <- readRDS("resources/var_ols_remember_author_cluster_101_tf.RDS")
redpill.var <- readRDS("resources/var_ols_remember_author_cluster_320_tf.RDS")
mental.var <- readRDS("resources/var_ols_remember_author_cluster_468_tf.RDS")
cod.var <- readRDS("resources/var_ols_remember_author_cluster_599_tf.RDS")
data.var <- readRDS("resources/var_ols_remember_author_cluster_631_tf.RDS")
realestate.var <- readRDS("resources/var_ols_remember_author_cluster_575_tf.RDS")
watches.var <- readRDS("resources/var_ols_remember_author_cluster_513_tf.RDS")
## toys.var <- readRDS("resources/var_ols_remember_author_cluster_23_tf.RDS")

## gun.var <- readRDS("resources/var_ols_remember_author_cluster_829_tf.RDS")

## art.var <- readRDS("resources/var_ols_remember_author_cluster_1636_tf.RDS")

## tanks.var <- readRDS("resources/var_ols_remember_author_cluster_1763_tf.RDS")

## tumblr.var <- readRDS("resources/var_ols_remember_author_cluster_1878_tf.RDS")

## energy.var <- readRDS("resources/var_ols_remember_author_cluster_1524_tf.RDS")

## diet.var <- readRDS("resources/var_ols_remember_author_cluster_1_tf.RDS")

if(!exists("breaks_pretty"))
  breaks_pretty <- pretty_breaks
options(xtable.floating = FALSE,
        xtable.timestamp = '',
        xtable.include.rownames=FALSE,
        math.style.negative=TRUE,
        booktabs = TRUE,
        xtable.format.args=list(big.mark=','),
        xtable.sanitize.text.function=identity
#        tikzDefaultEngine='xetex'
        )


# names(seattle.med.mu) <- unique(seattle.fit.df$subreddit) 
# idx <- which(names(seattle.med.mu)=='mariners')

## mariners.median.mu <- unlist(seattle.mu.stats[stat=='50%'][[paste0('mu[',idx,']')]])
## mariners.shock.size <- exp(mariners.median.mu + 1) - exp(mariners.median.mu)

## names(design.med.mu) <- unique(design.fit.df$subreddit)
## idx <- which(names(design.med.mu)=='onejob')
## onejob.median.mu <- unlist(design.mu.stats[stat=='50%'][[paste0('mu[',idx,']')]])
## onejob.shock.size <- exp(onejob.median.mu + 1) - exp(onejob.median.mu)

knit_hooks$set(document = function(x) {
    sub('\\usepackage[]{color}',
        '\\usepackage[]{color}', x, fixed = TRUE)
})
opts_chunk$set(fig.path="figures/knitr-")
opts_chunk$set(dev='pdf')
opts_chunk$set(external=TRUE)
opts_chunk$set(cache=FALSE)
overwrite <- FALSE

## shift_legend2 <- function(p) {
##   # check if p is a valid object
##   if(!(inherits(p, "gtable"))){
##     if(inherits(p, "ggplot")){
##       gp <- ggplotGrob(p) # convert to grob
##     } else {
##       message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
##       return(p)
##     }
##   } else {
##     gp <- p
##   }

##   # check for unfilled facet panels
##   facet.panels <- grep("^panel", gp[["layout"]][["name"]])
##   empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]), 
##                                USE.NAMES = F)
##   empty.facet.panels <- facet.panels[empty.facet.panels]

##   if(length(empty.facet.panels) == 0){
##     message("There are no unfilled facet panels to shift legend into. Returning original plot.")
##     return(p)
##   }

##   # establish name of empty panels
##   empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
##   names <- empty.facet.panels$name

##   # return repositioned legend
##   reposition_legend(p, 'center', panel=names)
## }

## shift_legend <- function(p){

##   # check if p is a valid object
##   if(!"gtable" %in% class(p)){

##   if("ggplot" %in% class(p)){
##       gp <- ggplotGrob(p) # convert to grob
##     } else {
##       message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
##       return(p)
##     }
##   } else {
##     gp <- p
##   }

##   # check for unfilled facet panels
##   facet.panels <- grep("^panel", gp[["layout"]][["name"]])
##   empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
##   empty.facet.panels <- facet.panels[empty.facet.panels]
##   if(length(empty.facet.panels) == 0){
##     message("There are no unfilled facet panels to shift legend into. Returning original plot.")
##     return(p)
##   }

##   # establish extent of unfilled facet panels (including any axis cells in between)
##   empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
##   empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
##                              max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
##   names(empty.facet.panels) <- c("t", "l", "b", "r")

##   # extract legend & copy over to location of unfilled facet panels
##   guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
##   if(length(guide.grob) == 0){
##     message("There is no legend present. Returning original plot.")
##     return(p)
##   }
##   gp <- gtable_add_grob(x = gp,
##                         grobs = gp[["grobs"]][[guide.grob]],
##                         t = empty.facet.panels[["t"]],
##                         l = empty.facet.panels[["l"]],
##                         b = empty.facet.panels[["b"]],
##                         r = empty.facet.panels[["r"]],
##                         name = "new-guide-box")

##   # squash the original guide box's row / column (whichever applicable)
##   # & empty its cell
##   guide.grob <- gp[["layout"]][guide.grob, ]
##   if(guide.grob[["l"]] == guide.grob[["r"]]){
##     gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
##   }
##   if(guide.grob[["t"]] == guide.grob[["b"]]){
##     gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
##   }
##   gp <- gtable_remove_grobs(gp, "guide-box")

##   return(gp)
## }

display.simulated.varmat <- function(){
    names <- sapply(1:dim(simulated_phi)[2],function(i) paste0('Y',i))
    obj = simulated_phi[,,]
    colnames(obj) <- names
    rownames(obj) <- names
    obj[obj < 0] <- paste0("\\color{yellow}",obj[obj < 0])
    obj[obj > 0] <- paste0("\\color{green}",obj[obj > 0])
    txt <- kable(obj,format='latex',booktabs=TRUE,escape=FALSE) #,row.names=TRUE,col.names=TRUE)
    return(txt)
}

display.recovered.varmat <- function(){

    stats <- simulation.phi.stats.95 
    median <- t(stats$medians[,,])
    upper <- t(stats$upper[,,])
    lower <- t(stats$lower[,,])
    names <- sapply(1:dim(median)[1],function(i) paste0('Y',i))

    sign.neg <- upper < 0
    sign.pos <- lower > 0
    obj <- round(median,2)
    obj[sign.neg] <- paste0("\\color{yellow}",obj[sign.neg])
    obj[sign.pos] <- paste0("\\color{green}",obj[sign.pos])
    colnames(obj) <- names
    rownames(obj) <- names
    
    txt <- kable(obj,format='latex',booktabs=TRUE,escape=FALSE) #,row.names=TRUE,col.names=TRUE)
    return(txt)
}

set.sig.level <- function(df,level){
    if(level=='95')
        old <- c('97.5%','2.5%')
    if(level=='90')
        old <- c('95%','5%')
    if(level=='85')
        old <- c('92.5%','7.5%')
    if(level=='80')
        old <- c('90%','10%')

    new <- c('upper','lower')
    
    if('data.table' %in% class(df)){
        setnames(df,old=old,new=new)
    } 
    
    if('list' %in% class(df)){
        names(df)[names(df)==old[1]] <- new[1]
        names(df)[names(df)==old[2]] <- new[2]
    }

    return(df)
}

#plot.irf.ols <- function(irf)

plot.irf <- function(df,level,ncol,level.b=NULL,level.c=NULL,level.names=NULL,override.level=NULL,scales="fixed",raw.irf=NULL,mu=NULL,include.all=FALSE){

    df <- copy(as.data.table(df))
    
    df[,":="(facet.title.tikz = paste0(name.col, " $\\rightarrow ", name.row),
             facet.title.expr = paste0(name.col, " %->% ", name.row))]

    if(opts_chunk$get("dev") == "tikz"){
        df[,facet.title := facet.title.tikz]
    } else {
        df[,facet.title := facet.title.expr]
    }

    keys <- c('level')

    if(is.null(raw.irf)){
        df.temp <- copy(df)
    } else {
       df.temp <- copy(raw.irf)
       df.temp[,":="(facet.title.tikz = paste0(name.col, " $\\rightarrow ", name.row),
                      facet.title.expr = paste0(name.col, " %->% ", name.row))]

        if(opts_chunk$get("dev") == "tikz"){
            df.temp[,facet.title := facet.title.tikz]
        } else {
            df.temp[,facet.title := facet.title.expr]
        }
    }
    

    if(!is.null(override.level)){
        df.temp <- set.sig.level(df.temp,override.level)
        if(include.all==FALSE){
            sig.pairs <- unique(df.temp[((upper < 0) |(lower >0)) & (name.row != name.col),.(facet.title)])
        } else {
            sig.pairs <- unique(df.temp[,.(facet.title)])
        }
    } else {
        set.sig.level(df.temp,level)
        if(include.all==FALSE){
            sig.pairs <- unique(df.temp[((upper < 0) |(lower >0)) & (name.row != name.col),.(facet.title)])
        } else {
            sig.pairs <- unique(df.temp[,.(facet.title)])
        }
    }

    if(!is.null(level.b)){
        df <- set.sig.level(df,level.b)
        setnames(df,old=c("lower","upper"),new=c("lower.b","upper.b"))
        keys <- c(keys,'level.b')
    }

    if(!is.null(level.c)){
        df <- set.sig.level(df,level.c)
        setnames(df,old=c("lower","upper"),new=c("lower.c","upper.c"))
        keys <- c(keys,'level.c')
    }

    df <- set.sig.level(df,level)

    if(!is.null(mu)){
        mu.df <- data.table(name.mu=unique(df$name.row),mu=mu)
    }

    df <- df[facet.title %in% sig.pairs$facet.title]

    if(!is.null(mu)){
        df <- df[mu.df,mu:=mu,on=.(name.row == name.mu)]
    }
        

    df[['median']] <- df[['50%']]
    p <- ggplot(df) #, aes(x=x,y=median,ymax=upper,ymin=lower))
    
    if(!is.null(level.c)){
        p <- p + geom_ribbon(aes(x=x,ymax=upper.c,ymin=lower.c,fill='level.c'))
        p <- p + geom_line(aes(x=x,y=upper.c,color='level.c'))
        p <- p + geom_line(aes(x=x,y=lower.c,color='level.c'))
    }
    if(!is.null(level.b)){
        p <- p + geom_ribbon(aes(x=x,ymax=upper.b,ymin=lower.b,fill='level.b'))
        p <- p + geom_line(aes(x=x,y=upper.b,color='level.b'))
        p <- p + geom_line(aes(x=x,y=lower.b,color='level.b'))
    }

    colors <- c('grey80','grey65','grey45')
    colors <- colors[1:length(keys)]
    names(colors) <- rev(keys)

    if(is.null(level.names))
        level.names <- keys

    p <- p + geom_ribbon(aes(x=x,ymax=upper,ymin=lower,fill='level'))
    p <- p + geom_line(aes(x=x,y=upper,color='level'))
    p <- p + geom_line(aes(x=x,y=lower,color='level'))

    p <- p + geom_line(aes(x=x,y=median),show.legend=FALSE)

    if(!is.null(mu)){
        p <- p + geom_hline(aes(yintercept=mu,linetype='dashed'),data=df)
        p <- p + scale_linetype_manual(name="",values=c("dashed"),breaks=c("dashed"),labels=("baseline"))
    }

    if(opts_chunk$get("dev") != "tikz"){
        p <- p + facet_wrap(.~facet.title, labeller=label_parsed, ncol=ncol,scales=scales)
    } else {
        p <- p + facet_wrap(.~facet.title, ncol=ncol,scales=scales)
    }
    p <- p + scale_fill_manual(name='Credible interval',values=colors,breaks=rev(keys),labels=rev(level.names))
    p <- p + scale_color_manual(name='Credible interval',values=colors,breaks=rev(keys),labels=rev(level.names))
    p <- p + xlab("Forecast week") + ylab("Impulse response")
    p <- p + theme(legend.position='bottom')
    p <- p + guides(fill=FALSE)
    return(p)

}

plot.coefs <- function(stats, names, stats.b=NULL, stats.c=NULL,stat.names=NULL,point.shape.name=NULL,include.all=F,ncols=1){
    median <- t(stats$medians[,,])
    upper <- t(stats$upper[,,])
    lower <- t(stats$lower[,,])

    if(!is.null(stats.b)){
        upper.b = t(stats.b$upper[,,])
        lower.b = t(stats.b$lower[,,])
    }

    if(!is.null(stats.c)){
        upper.c = t(stats.c$upper[,,])
        lower.c = t(stats.c$lower[,,])
    }

    plot.df <- list()
    m <- dim(stats$median[,,])[1]#$
    idx <- 1
    for(i in 1:m){
        for(j in 1:m){
            row <- list('i'=names[[i]],'j'=names[[j]],'coef'=median[i,j],'upper'=upper[i,j],'lower'=lower[i,j],name.expr=paste0(names[[i]],"%->%",names[[j]]),name.tex=paste0(names[[i]],"$\\rightarrow$",names[[j]]))

            if(!is.null(stats.b)){
                row['upper.b'] = upper.b[i,j]
                row['lower.b'] = lower.b[i,j]
            }

            if(!is.null(stats.c)){
                row['upper.c'] = upper.c[i,j]
                row['lower.c'] = lower.c[i,j]
            }

            plot.df[[idx]] <- row
            idx <- idx + 1
        }
    }

    plot.df <- rbindlist(plot.df)

    if(include.all == F){
        plot.df <- plot.df[((upper < 0) | (lower > 0)) & (i!=j)]
    }

    plot.df <- plot.df[,':='(name.tex = factor(name.tex,levels=rev(name.tex)))]
    plot.df <- plot.df[,':='(name.expr = factor(name.expr,levels=rev(name.expr)))]

    if(opts_chunk$get("dev") == 'tikz'){
        plot.df$name = plot.df$name.tex
    } else {
        plot.df$name = plot.df$name.expr
    }

    if(ncols > 1){
        groups <- unique(plot.df[,.(name)])
        groups <- groups[,group := cut(1:nrow(groups),ncols)]
        plot.df <- plot.df[groups,on=.(name=name)]
    }

    p <- ggplot(plot.df, aes(x=name,y=coef,ymin=lower,ymax=upper))

    keys <- c()

    if(!is.null(stats.c)){
        p <- p + geom_segment(data=plot.df,aes(x=name,xend=name,y=lower.c,yend=upper.c,color='stats.c'),size=0.8)
        keys <- c(keys,'stats.c')
    }


    if(!is.null(stats.b)){
        p <- p + geom_segment(data=plot.df,aes(x=name,xend=name,y=lower.b,yend=upper.b,color='stats.b'),size=0.8)
        keys <- c(keys,'stats.b')
    }

    p <- p + geom_segment(aes(x=name,xend=name,y=lower,yend=upper,color='stats'),size=0.8)
    keys <- c(keys,'stats')

    if (!is.null(point.shape.name)){
        p <- p + geom_point(aes(shape=point.shape.name),size=1.5)
    } else {
        p <- p + geom_point(size=1.5)
    }

    if(ncols > 1){
        p <- p + facet_wrap(group~.,scales='free_y',nrow=1)
    }

    p <- p + coord_flip()

    if(opts_chunk$get("dev") == 'tikz'){
        p <- p + scale_x_discrete(labels = levels(plot.df$name),plot.df$name)#$
    } else {
        p <- p + scale_x_discrete(labels = parse(text=levels(plot.df$name)),breaks=levels(plot.df$name))#$
    }


    p <- p + ylab("") + xlab("Coefficient")
#    p <- p + scale_size_manual(name='Credible interval',values=c(stats=1,stats.b=1.25,stats.c=1.5),breaks=c("stats","stats.b","stats.c"),labels=stat.names)

    colors <- c('grey45','grey65','grey80')
    colors <- colors[1:length(keys)]
    names(colors) <- rev(keys)

    if(is.null(stat.names))
        stat.names <- keys

    p <- p + scale_color_manual(name='Credible interval',values=colors,breaks=rev(keys),labels=stat.names)
    if(ncols >1){
        p <- p + theme(legend.position='bottom',panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),strip.background=element_blank(),strip.text.x=element_blank())
    } else {
        p <- p + theme(legend.position='bottom',panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())
    }
    return(p)
}

plot.irf.ols.data <- function(irf){
    irf <- irf[,":="(impulse.subreddit = paste0("r.",gsub("^X","",impulse.subreddit)),
                     response.subreddit = paste0("r.",gsub("^X","",response.subreddit)))]
    
    irf <- irf[,":="(facet.title.tikz = paste0(impulse.subreddit, " $\\rightarrow ", response.subreddit),
                     facet.title.expr = paste0(impulse.subreddit, " %->% ", response.subreddit))]

    sig.pairs <- unique(irf[((irf.upper < 0) |(irf.lower >0)) & (impulse.subreddit != response.subreddit),.(facet.title.expr)])

    irf <- irf[facet.title.expr %in% sig.pairs$facet.title.expr]

    return(irf)
}

plot.ols.forecast.data <- function(var.ypred,ar.ypred){
    var.ypred[,model:='var']
    ar.ypred[,model:='ar']

    plot.data <- rbind(var.ypred,ar.ypred)
    plot.data[,subreddit:=gsub("^X","",subreddit)]

    return(plot.data)
}

plot.ols.forecast <- function(var.ypred, ar.ypred){
    plot.data <- plot.ols.forecast.data(var.ypred, ar.ypred)
    p <- ggplot(plot.data,aes(x=week,y=y))
    p <- p + geom_line(color='black')
    p <- p + geom_line(aes(x=week,y=y.fit),color='red',linetype='dashed')
    p <- p + geom_line(aes(x=week,y=y.pred),color='red',linetype='dashed')
    p <- p + geom_ribbon(aes(x=week,ymax=y.pred.upper,ymin=y.pred.lower),alpha=0.8)
    p <- p + facet_grid(subreddit ~ model)
    return(p)
}

plot.coef.ols.data <- function(coef){

    proc.coefs <- function(coef,i){
        tab <- coef[[i]]
        nms <- rownames(tab)
        tab <- data.table(tab)
        tab <- tab[,term:=nms]
        commenseterms <- tab$term[grepl(".*\\.l1",tab$term)]
        coef <- tab[,"Estimate",with=F]
        std.err <- tab[,"Std. Error",with=F]

        ci.upper <- coef + 1.98 * std.err
        ci.lower <- coef - 1.98 * std.err
        tab[,ci.upper := ci.upper]
        tab[,ci.lower := ci.lower]
        tab[term %in% commenseterms,j:=gsub("\\.l1","",gsub("^X","",term))]
        tab[,i:=gsub("^X","",i)]
        tab[,name.tikz:= paste0(j, " $\\rightarrow$ ", i)]
        tab[,name.expr:= paste0("r.", j, " %->% ", "r.", i)]
        return(tab)
    }

    tab <- rbindlist(lapply(names(coef),function(i) proc.coefs(coef,i)))
    ctab <- tab[!is.na(j)]
    ctab[,model:=1]

    ctab[,estimate:=Estimate]
    ctab[["std.error"]] <- ctab[["Std. Error"]]
    setnames(ctab,old=c('term'),new=c('og.term'))
    setnames(ctab,old=c('name.expr'),new=c('term'))
                                        #ctab[,term:=sapply(term,function(c) str2expression(c))]
    ctab[,term:=as.factor(term)]
    ctab <- ctab[i!=j]

    return(ctab)
}


plot.simulated.forecast <- function(level='95'){
    level='95'

    simulation.fit.df <- as.data.table(simulation.fit.df)[,type:='fit']
    simulation.forecast.df <- simulation.forecast.df[,type:='forecast']
    df <- rbind(simulation.fit.df,simulation.forecast.df)

    if(opts_chunk$get("dev") == 'tikz')
        df <- df[,subreddit := paste0(gsub('y','$y_',subreddit),'$')]

    simulation.forecast.data <- set.sig.level(simulation.forecast.data,level)
    
    if(opts_chunk$get("dev") == 'tikz')
        simulation.forecast.data <- simulation.forecast.data[,subreddit := paste0(gsub('y','$y_',subreddit),'$')]

    simulation.forecast.data <- simulation.forecast.data[,median:=simulation.forecast.data[['50%']]]
    max.x <- max(df[type=='fit',x])
    p <- ggplot()
    p <- p + geom_line(aes(x=week,y=N_authors),data=df,color='black')
    p <- p + geom_line(aes(x=week,y=mean),data=simulation.forecast.data[type=='forecast'],color='red',linetype='dashed')
    p <- p + geom_ribbon(aes(x=week,ymax=upper,ymin=lower),data=simulation.forecast.data,alpha=0.4)
#    p <- p + scale_linetype_manual(limits=c("solid","dashed"),values=c("fit","forecast"),breaks=c("fit",'forecast'))
    p <- p + facet_wrap(.~subreddit,scales='free_y',nrow=2,strip.position='left')
    p <- p + theme(legend.position='left')
    p <- p + ylab('') + xlab('week')
    p

}

plot.wallpaper.ifs <- function(level='80'){

    df <- copy(as.data.table(wallpaper.irf.ortho.data))
    
    df[,":="(facet.title.tikz = paste0(name.row, " $\\rightarrow ", name.col),
             facet.title.expr = paste0(name.row, " %->% ", name.col))]

    if(opts_chunk$get("dev") == "tikz"){
        df[,facet.title := facet.title.tikz]
    } else {
        df[,facet.title := facet.title.expr]
    }

    df <- set.sig.level(df,level)

    sig.pairs <- unique(df[((upper < 0) |(lower >0)) & (name.row != name.col),.(facet.title)])

    df <- df[facet.title %in% sig.pairs$facet.title]
        
    df[['median']] <- df[['50%']]

    p <- ggplot(df, aes(x=x,y=median,ymax=upper,ymin=lower)) + geom_line() + geom_ribbon(alpha=0.5)

    ## my_labeller <- as_labeller(function(value){
    ##     bquote(.(parse(text=value)))
    ## })

    p <- p + facet_wrap(.~facet.title,labeller=label_parsed)
    p <- p + xlab("Forecast week") + ylab("Impulse response")
    p


}


plot.simulation.irfs <- function(){
    
    df <- copy(as.data.table(simulation.irf.forecast.data))
    
    if(opts_chunk$get("dev") == 'tikz'){
        df <- df[,name.row := paste0(gsub('y','$y_',name.row),'$')]
        df <- df[,name.col := paste0(gsub('y','$y_',name.col),'$')]
        focal.sub <- '$y_1$'

    } else {
        focal.sub <- 'y1'
    }

    matnames <- unique(df$name.row)

    level <- '95'

    if(level=='95')
        setnames(df,old=c('97.5%','2.5%'),new=c('upper','lower'),skip_absent=T)
    if(level=='90')
        setnames(df,old=c('95%','5%'),new=c('upper','lower'),skip_absent=T)
    if(level=='80')
        setnames(df,old=c('85%','15%'),new=c('upper','lower'),skip_absent=T)

    df <- df[(name.row == focal.sub) | (name.col==focal.sub)]
    df <- df[(name.row != name.col)]

    if(opts_chunk$get("dev") == 'tikz'){
        df <- df[,":="(name.row = factor(name.row,levels=matnames),
                                         name.col = factor(name.col,levels=matnames),
                                        dir=ifelse(name.row==focal.sub,paste0("$y_m\\rightarrow",focal.sub,'$'),paste0("$",focal.sub,"\\rightarrow","y_m$")),
                                         name=ifelse(name.row==focal.sub,name.col,name.row)
                                     )]
    } else {
        df <- df[,":="(name.row = factor(name.row,levels=matnames),
                        name.col = factor(name.col,levels=matnames),
                        dir=ifelse(name.row==focal.sub,paste0("on ",focal.sub),paste0("of ",focal.sub)),
                        name=ifelse(name.row==focal.sub,name.col,name.row)
                )]
    }

    df[['median']] <- df[['50%']]

    p <- ggplot(df, aes(x=x, y=median, ymin=lower,ymax=upper)) + geom_line() + geom_ribbon(alpha=0.5) + facet_grid(dir~name)
    
    p <- p + xlim(1,8)
    p <- p + xlab("Forecast week") + ylab("Impulse response")
    p

}

plot.coef.barchart <- function(plot.df){
    plot.df <- plot.df[,':='(name = factor(name,levels=rev(name)))]
    p <- ggplot(plot.df, aes(x=name,y=coef,ymin=lower,ymax=upper,color='estimate',shape='estimate')) + geom_pointrange() + geom_point(aes(x=name,y=true,color='value',shape='value'),size=3) + coord_flip()
    p <- p + scale_color_manual(name='',values=c('black','blue'),breaks=c("estimate","value"),labels=c("Estimate","True value"))
    p <- p + scale_shape_manual(name='',values=c(estimate=19,value=1), breaks=c("estimate","value"), labels=c("Estimate","True value"))
    p <- p + ylab("Value") + xlab("Coefficient") + theme(legend.position='top')

}

recovered.params.barchart <- function(){

    stats <- simulation.phi.stats.95 
    median <- t(stats$medians[,,])
    upper <- t(stats$upper[,,])
    lower <- t(stats$lower[,,])
    names <- sapply(1:dim(median)[1],function(i) paste0('Y',i))

    plot.df <- list()
    m <- dim(stats$median[,,])[1]#$
    idx <- 1
    for(i in 1:m){
        for(j in 1:m){
            plot.df[[idx]] <- list('i'=i,'j'=j,'coef'=median[i,j],'upper'=upper[i,j],'lower'=lower[i,j],'true'=simulated_phi[,,][i,j])
            idx <- idx + 1
        }
    }
    plot.df <- rbindlist(plot.df)
    if(opts_chunk$get("dev") == "tikz"){#$
        plot.df <- plot.df[,name:=paste0("\\phi_{",i,',',j,'}')]
    } else {
        plot.df <- plot.df[,name:=paste0("phi[list(",i,',',j,')]')]
    }
   
    p <- plot.coef.barchart(plot.df)
    
    if(opts_chunk$get("dev") == 'tikz'){
    p <- p + scale_x_discrete(labels = plot.df$name)#$
    } else {
    p <- p + scale_x_discrete(labels = parse(text=levels(plot.df$name)))#$
    }
    print(p)

}

wallpaper.barchart <- function(){
    
    stats <- wallpaper.phi.stats.80
    median <- t(stats$medians[,,])
    upper <- t(stats$upper[,,])
    lower <- t(stats$lower[,,])
    
    names <- unique(wallpaper.fit.df$subreddit)
    plot.df <- list()
    m <- dim(stats$median[,,])[1]#$
    idx <- 1
    for(i in 1:m){
        for(j in 1:m){
            plot.df[[idx]] <- list('i'=names[[i]],'j'=names[[j]],'coef'=median[i,j],'upper'=upper[i,j],'lower'=lower[i,j],name.expr=paste0(names[[i]],"%->%",names[[j]]),name.tex=paste0(names[[i]],"$\\rightarrow$",names[[j]]))
            idx <- idx + 1
        }
    }

    plot.df <- rbindlist(plot.df)
    plot.df <- plot.df[,':='(name.tex = factor(name.tex,levels=rev(name.tex)))]
    plot.df <- plot.df[,':='(name.expr = factor(name.expr,levels=rev(name.expr)))]

    if(opts_chunk$get("dev") == 'tikz'){
        plot.df$name = plot.df$name.tex
    } else {
        plot.df$name = plot.df$name.expr
    }

    plot.df <- plot.df[((upper < 0) | (lower > 0)) & (i!=j)]

    p <- ggplot(plot.df, aes(x=name,y=coef,ymin=lower,ymax=upper,color='estimate',shape='estimate')) + geom_pointrange() + coord_flip()
    p <- p + ylab("Value") + xlab("Coefficient") + theme(legend.position='top')


    if(opts_chunk$get("dev") == 'tikz'){
        p <- p + scale_x_discrete(labels = as.character(plot.df$name))#$
    } else {
        p <- p + scale_x_discrete(labels = parse(text=as.character(plot.df$name)))#$
    }

    print(p)
}

display.estimated.varmat <- function(){

    stats <- wallpaper.phi.stats.80
    
    median <- t(stats$medians[,,])
    upper <- t(stats$upper[,,])
    lower <- t(stats$lower[,,])
    sign.neg <- upper < 0
    sign.pos <- lower > 0
    obj <- round(median,2)
    names <- unique(wallpaper.fit.df$subreddits)
    obj[sign.neg] <- paste0("\\color{yellow}",obj[sign.neg])
    obj[sign.pos] <- paste0("\\color{green}",obj[sign.pos])
    colnames(obj) <- names
    rownames(obj) <- names
    
    txt <- kable(obj,format='latex',booktabs=TRUE,escape=FALSE) #,row.names=TRUE,col.names=TRUE)
    return(txt)
}

#    p <- plot.ortho.irf(design.irf.ortho.units.data[x<=20],level='80',ncol=2,'90','95',c("80%","90%","95%"),"95",scales='free_y',raw.irf = design.irf.ortho.data[x<=20],mu=design.med.mu)
#    names <- unique(design.fit.df$subreddit)
#    df <- data.table(name=names,mu=exp(design.med.mu))
#    p <- p + geom_hline(aes(yintercept=mu),data=df,linetype='dashed',color='grey10')
#    p <- p + guides(color=guide_legend(ncol=1,byrow=TRUE))
#    p <- p + theme(legend.position='bottom',legend.direction='horizontal')
##    sr.names <- unique(design.irf.forecast.units.data[,.(name.row)])
##     sr.names <- sr.names[,part := cut(1:nrow(sr.names),2)]
##     parts <- unique(sr.names$part)
    
## p <- plot.irf(design.irf.forecast.units.data[x<=20][name.row%in%[sr.names[part==parts[[1]],name.row]]],'80',ncol=6,'90','95',c("80%","90%","95%"),"95",scales='free_y',raw.irf = design.irf.forecast.data[x<=20][name.row%in%parts[[1]]],mu=design.med.mu,include.all=TRUE)
##    sr.names <- unique(seattle.irf.forecast.units.data[,.(name.col)])
##     sr.names <- sr.names[,part := cut(1:nrow(sr.names),3)]
##     sr.names <- sr.names[,mu:=seattle.med.mu]
##     parts <- unique(sr.names$part)
    
##     p <- plot.irf(seattle.irf.forecast.units.data[x<=20][name.col%in%sr.names[part==parts[[1]],name.col]],'80',ncol=6,'90','95',c("80%","90%","95%"),"95",scales='free_y',raw.irf =seattle.irf.forecast.data[x<=20][name.col%in%sr.names[part==parts[[1]],name.col]],mu=seattle.med.mu,include.all=TRUE)
##     p <- p + guides(color=guide_legend(ncol=2,byrow=TRUE,order=1),linetype=guide_legend(),order=2)
## p <- p + theme(legend.position='bottom',legend.box='horizontal')
## p

 ## stats <- copy(simulation.phi.stats.95)
 ##    median <- t(stats$medians[,,])
 ##    upper <- t(stats$upper[,,])
 ##    lower <- t(stats$lower[,,])
 ##    names <- sapply(1:dim(simulated_phi)[2],function(i) paste0('Y',i))
 ##    plot.df <- list()
 ##    m <- dim(stats$median[,,])[1]#$
 ##    idx <- 1
 ##    for(i in 1:m){
 ##        for(j in 1:m){
 ##            plot.df[[idx]] <- list('i'=i,'j'=j,'coef'=median[i,j],'upper'=upper[i,j],'lower'=lower[i,j],
 ##            'true'=simulated_phi[,,][i,j],
 ##            name.expr=paste0(names[[i]],"%->%",names[[j]]),
 ##            name.tex=paste0(names[[i]],"$\\rightarrow$",names[[j]]))
 ##            idx <- idx + 1
 ##        }
 ##    }
 ##    plot.df <- rbindlist(plot.df)

 ##    plot.df <- plot.df[,':='(name.tex = factor(name.tex,levels=rev(name.tex)))]
 ##    plot.df <- plot.df[,':='(name.expr = factor(name.expr,levels=rev(name.expr)))]

 ##    if(opts_chunk$get("dev") == 'tikz'){
 ##        plot.df$name = plot.df$name.tex
 ##    } else {
 ##        plot.df$name = plot.df$name.expr
 ##    }

 ##    N.sim.wrong <- nrow(plot.df[(true < lower) | (true > upper)])
    
 ##    p <- plot.coefs(simulation.phi.stats.80, names, simulation.phi.stats.90, simulation.phi.stats.95, stat.names=c("80%","90%","95%"),include.all=T,point.shape.name='estimate',ncols=2)

 ##    p <-  p + geom_point(aes(x=name,y=true,shape='value'),data=plot.df,size=1.5) + coord_flip()


 ##    p <- p + scale_shape_manual(name='',values=c(estimate=19,value=1), breaks=c("estimate","value"), labels=c("Estimate","True value"))
 ##    p <- p + ylab("Value") + xlab("Coefficient") + theme(legend.position='bottom',panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
 ##    p <- p + scale_x_discrete(labels = parse(text=levels(plot.df$name)),breaks=plot.df$name)#$
 ##    p <- p + guides(color=guide_legend(direction='vertical',ncol=2,byrow=F), shape=guide_legend(direction='vertical'))
 ##    print(p)
