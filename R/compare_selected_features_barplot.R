### selected features in each patients
features = c("3702")
library(gridExtra)
compare_selected_features_barplot <- function(features,dds.lym,lym.norm,ecz.norm,dds.ecz){
    mat.norm <- lym.norm
    dds <- dds.lym
    if (length(features) == 1){
        my.features <- data.frame(cbind(mat.norm[which(rownames(mat.norm)%in%features), ],
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    } else {
    my.features <- data.frame(cbind(t(mat.norm[which(rownames(mat.norm)%in%features), ]),
                                    as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                    as.matrix(colData(dds)$ID)))
    }
    names(my.features) <- c(features, "condition", "pedigree", "ID")
    my.features <- my.features[order(my.features$pedigree, my.features$condition), ]
    plot.mat_lym <- melt(my.features, id=c("condition", "pedigree", "ID"))
    plot.mat_lym$value <- as.numeric(plot.mat_lym$value)
    plot.mat_lym$class <- rep("lym", dim(plot.mat_lym)[1])
    plot.mat_lym$condition <- paste0(plot.mat_lym$class ,plot.mat_lym$condition)
    plot.mat_lym$class <- NULL
    #-------------------------------------------------------------------------------
    mat.norm <- ecz.norm
    dds <- dds.ecz
    if (length(features) == 1){
        my.features <- data.frame(cbind(mat.norm[which(rownames(mat.norm)%in%features), ],
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    }else{
        my.features <- data.frame(cbind(t(mat.norm[which(rownames(mat.norm)%in%features), ]),
                                    as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                    as.matrix(colData(dds)$ID)))
    }
    names(my.features) <- c(features, "condition", "pedigree", "ID")
    my.features <- my.features[order(my.features$pedigree, my.features$condition), ]
    plot.mat_ecz <- melt(my.features, id=c("condition", "pedigree", "ID"))
    plot.mat_ecz$value <- as.numeric(plot.mat_ecz$value)
    plot.mat_ecz$class <- rep("ecz", dim(plot.mat_ecz)[1])
    plot.mat_ecz$condition <- paste0(plot.mat_ecz$class ,plot.mat_ecz$condition)
    plot.mat_ecz$class <- NULL
#-------------------------------------------------------------------------------
    plot_mat <- rbind(plot.mat_lym,plot.mat_ecz)
    min_lym <- plot_mat %>% filter(plot_mat$condition %in% c("lymI","lymnI")) %>%
        group_by(pedigree,variable) %>% summarise(values = sum(value)) %>%
        group_by(variable) %>% summarize(min(values))

    plot_mat_ratio_lym <- plot.mat_lym %>%
        group_by(pedigree) %>%
        arrange(condition) %>%
        summarise(value = value[condition == "lymI"] / value[condition == "lymnI"] )
    plot_mat_ratio_lym$class <- "lym"
    plot_mat_ratio_ecz <- plot.mat_ecz %>%
        group_by(pedigree) %>%
        arrange(condition) %>%
        summarise(value = value[condition == "eczI"] / value[condition == "ecznI"] )
    plot_mat_ratio_ecz$class <- "ecz"

    plot_mat_ratio <- rbind(plot_mat_ratio_lym,plot_mat_ratio_ecz)
    plot_mat$condition <- factor(plot_mat$condition, levels = c("lymI", "lymnI","eczI", "ecznI" ))
    bar.cols <- brewer.pal(nlevels(factor(plot_mat$condition)), "Set1")
    p1 <- ggplot(plot_mat, aes(pedigree, value, fill=condition)) +
        geom_col() +
        geom_hline(yintercept = min_lym$`min(values)`)+
        facet_grid(variable~.)+
        scale_fill_manual(values=c(bar.cols))+
        labs(title="") + xlab("") + ylab("") +
        theme_bw() + theme(legend.title=element_blank())+
        theme(legend.position = c(1, 1), legend.justification = c(1, 1),
              legend.direction = "horizontal", legend.key.size = unit(.3,units = "cm"),
              legend.background = element_rect(color = "white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())

    p2 <- ggplot(plot_mat_ratio, aes(pedigree, value, fill=class)) +
        geom_col() +
        scale_fill_manual(values=c(bar.cols))+
        labs(title="lesion - nonlesion ratio") + xlab("") + ylab("") +
        theme_bw() + theme(legend.title=element_blank())+
        theme(legend.position = c(1, 1), legend.justification = c(1, 1),
              legend.direction = "horizontal", legend.key.size = unit(.3,units = "cm"),
              legend.background = element_rect(color = "white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())



    mat.norm <- lym.norm
    dds <- dds.lym
    if (length(features) == 1){
        my.features <- data.frame(cbind(scale(mat.norm[which(rownames(mat.norm)%in%features), ]),
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    } else {
        my.features <- data.frame(cbind(t(scale(mat.norm[which(rownames(mat.norm)%in%features), ])),
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    }
    names(my.features) <- c(features, "condition", "pedigree", "ID")
    my.features <- my.features[order(my.features$pedigree, my.features$condition), ]
    plot.mat_lym <- melt(my.features, id=c("condition", "pedigree", "ID"))
    plot.mat_lym$value <- as.numeric(plot.mat_lym$value)
    plot.mat_lym$class <- rep("lym", dim(plot.mat_lym)[1])
    plot.mat_lym$condition <- paste0(plot.mat_lym$class ,plot.mat_lym$condition)
    plot.mat_lym$class <- NULL
    #-------------------------------------------------------------------------------
    mat.norm <- ecz.norm
    dds <- dds.ecz
    if (length(features) == 1){
        my.features <- data.frame(cbind(scale(mat.norm[which(rownames(mat.norm)%in%features), ]),
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    }else{
        my.features <- data.frame(cbind(t(scale(mat.norm[which(rownames(mat.norm)%in%features), ])),
                                        as.matrix(colData(dds)$condition), as.matrix(colData(dds)$pedigree),
                                        as.matrix(colData(dds)$ID)))
    }
    names(my.features) <- c(features, "condition", "pedigree", "ID")
    my.features <- my.features[order(my.features$pedigree, my.features$condition), ]
    plot.mat_ecz <- melt(my.features, id=c("condition", "pedigree", "ID"))
    plot.mat_ecz$value <- as.numeric(plot.mat_ecz$value)
    plot.mat_ecz$class <- rep("ecz", dim(plot.mat_ecz)[1])
    plot.mat_ecz$condition <- paste0(plot.mat_ecz$class ,plot.mat_ecz$condition)
    plot.mat_ecz$class <- NULL
    #-------------------------------------------------------------------------------
    plot_mat <- rbind(plot.mat_lym,plot.mat_ecz)
    min_lym <- plot_mat %>% filter(plot_mat$condition %in% c("lymI","lymnI")) %>%
        group_by(pedigree,variable) %>% summarise(values = sum(value)) %>%
        group_by(variable) %>% summarize(min(values))

    plot_mat_ratio_lym <- plot.mat_lym %>%
        group_by(pedigree) %>%
        arrange(condition) %>%
        summarise(value = value[condition == "lymI"] / value[condition == "lymnI"] )
    plot_mat_ratio_lym$class <- "lym"
    plot_mat_ratio_ecz <- plot.mat_ecz %>%
        group_by(pedigree) %>%
        arrange(condition) %>%
        summarise(value = value[condition == "eczI"] / value[condition == "ecznI"] )
    plot_mat_ratio_ecz$class <- "ecz"

    plot_mat_ratio <- rbind(plot_mat_ratio_lym,plot_mat_ratio_ecz)
    plot_mat$condition <- factor(plot_mat$condition, levels = c("lymI", "lymnI","eczI", "ecznI" ))
    bar.cols <- brewer.pal(nlevels(factor(plot_mat$condition)), "Set1")
    p3 <- ggplot(plot_mat, aes(pedigree, value, fill=condition)) +
        geom_col() +
        geom_hline(yintercept = min_lym$`min(values)`)+
        facet_grid(variable~.)+
        scale_fill_manual(values=c(bar.cols))+
        labs(title="scaling lym and ecz separately") + xlab("") + ylab("") +
        theme_bw() + theme(legend.title=element_blank())+
        theme(legend.position = c(1, 1), legend.justification = c(1, 1),
              legend.direction = "horizontal", legend.key.size = unit(.3,units = "cm"),
              legend.background = element_rect(color = "white"))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())

    grid.arrange(p1, p2,p3, nrow = 3)
}

