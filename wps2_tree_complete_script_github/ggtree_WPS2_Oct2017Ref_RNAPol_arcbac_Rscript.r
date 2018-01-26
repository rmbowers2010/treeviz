#### #### 	GGTREE SCRIPT TO PLOT CONSERVED MARKER TREE WITH ASSOCIATED METADATA 	#### ####

#!/usr/bin/Rscript


# this line is needed to run Rscript from command line
args <- commandArgs(trailingOnly = TRUE)


library(ggtree)
library(ape)
library(phylotools)
library(grid)
library(geiger) # tips() function
library(phytools)
library(tidyr)
library(phyloseq)

#	args[1] = tree_file.nwk
#	args[2] = mapping_file.txt
#	args[3] = circular or linear
#	args[4] = circular or linear

#	args[ ] = lookup.tab
#	args[ ] = reference_lookup.tab
#	args[ ] = cog/function table, maybe??
#	args[ ] = taxa_to_color.txt

# one issue -> figure out how to more efficiently load metadata
# add query map file, then code the ref as 'ref' in 'sample_type' column

 	tree_file <- read.tree(args[1])
	#tree_file <- read.tree(list.files(pattern="FTWAG"))

# sample mapping file
# CURRENTLY: setup to have first two columns be the important columns
#	Col1: TAXON_OID, Col2: label you want on tree, needs to be unique
	mapping_file <- read.table(args[2], sep="\t", header=TRUE)

	# adding concat.faa.faa to mapping_file table
	tmp <- c("concatFaa", "WPS-2_concatenated", NA, NA, NA, NA, NA, NA, NA, NA, NA)
	tmp.df <- t(data.frame(tmp))
	colnames(tmp.df) <- colnames(mapping_file)
	tmp.df <- as.data.frame(tmp.df)
	rownames(tmp.df) <- NULL
	mapping_file <- rbind(mapping_file, tmp.df)	



	#mapping_file <- read.table("../../wps-2_taxids_and_metadata_rmDups.txt", sep="\t", header=TRUE)
	mapping_file <- mapping_file[,c(1,2)]; colnames(mapping_file) <- c("V1", "V2")
	mapping_file$V3 <- c("WPS2")

# reference lookup
	reference_lookup <- read.table("/projectb/scratch/rmbowers//DBASES_USEFUL_TO_NEW_LIFE_PROJECTS/arcbac27October2017_genomes_fna_faa/bacteria27October2017.lookup", header=FALSE, sep="\t")
	reference_lookup$V1 <- gsub("IMG", "", reference_lookup$V1)
	reference_lookup$V3 <- c("ref")

# tree_map 
	tree_map <- rbind(mapping_file, reference_lookup)
	tree_map <- tree_map[tree_map$V1 %in% tree_file$tip.label, ]
	tree_map <- tree_map[!duplicated(tree_map$V1),]
	tree_map$V2 <- paste(tree_map$V2, tree_map$V1, sep = " | ")
	tree_map$Samples_added_to_Tree <- 1; tree_map$Samples_added_to_Tree [tree_map$V3 =="ref"] <- 0


# Create otu_tab from mapping_file
	otu_tab <- tree_map %>% select(V2, Samples_added_to_Tree)
	rownames(otu_tab) <- otu_tab[,1]; otu_tab[,1] <- NULL
	OTU <- otu_table(as.matrix(otu_tab), taxa_are_rows=TRUE)


# swap tree_file tip labels - these tips use 'sample_lookup' from genomestoreferencetree
# depending on version of phylotools: sub.taxa.label might be needed instead of sub.tip.label
# looks like some of the bad genomes were ditched, ex: Q0DHDNLB=AHXAA=0%complete, 
	tree_file <- sub.taxa.label(tree_file, tree_map[,c("V1", "V2")])

# outgroup
	outgroup_names <- tree_file$tip.label[grepl(paste(c("Parcu","Microgenomates"), collapse="|"), tree_file$tip.label)]
	outgroup_fullnames <- tree_map[tree_map$V2 %in% outgroup_names, ]$V2
	root_nodeID <- MRCA(tree_file, outgroup_fullnames)
	tree_rerooted <- reroot(tree_file, root_nodeID)
	root_nodeID_rerootID <- MRCA(tree_rerooted, outgroup_fullnames)
	tree_rerooted_rerooted <- reroot(tree_rerooted, root_nodeID_rerootID)
	tree_rerooted_rerooted <- rotate(tree_rerooted_rerooted, root_nodeID_rerootID)



# drop tips by long branch length if needed
	# http://blog.phytools.org/2012/08/trick-to-name-edges-of-phylo-object-for.html
	# IF THERE ARE NODES WITH TWO LONG BRANCHES, CHANCES ARE THE SMALLER WILL NOT HAVE A LONG BRANCH UNTIL PLOTTING ON ITS OWN, MEANING, YOU MIGHT HAVE TO REPEAT THIS ITERATIVELY
	X <- tree_rerooted_rerooted$edge
	# replaces second column in edge table with tip labels
	X[X[,2]%in%1:length(tree_rerooted_rerooted$tip),2]<- tree_rerooted_rerooted$tip[X[X[,2]%in%1:length(tree_rerooted_rerooted$tip),2]]
	# name edge.lengths with corresponding tip labels based on the indexing done in X
	#names(tree_rerooted_rerooted$edge.length)<-paste(X[,1],X[,2],sep=",")
	names(tree_rerooted_rerooted$edge.length)<-X[,2]
	tipNames_to_drop <- names(head(tree_rerooted_rerooted$edge.length[order(-tree_rerooted_rerooted$edge.length)], 50))
	tree_rerooted_rerooted <- drop.tip(tree_rerooted_rerooted, tipNames_to_drop)

# drop tips by name if you need drop any, using drop.tip() function

	tips_to_drop_manual <- c("2504756006","2534681894","2548877146","259333924","2606217800","2695420353","2547132491")
	tips_to_drop_manual <- tree_rerooted_rerooted$tip.label[grepl(paste(tips_to_drop_manual, collapse="|"), tree_rerooted_rerooted$tip.label)]


	tree_rerooted_rerooted <- drop.tip(tree_rerooted_rerooted, tips_to_drop_manual)


# COLORING IS OFF, MISSING FIRST ORGANISM ("ACIDOBACTERIA") AND MCRA IS PULLING THE WRONG NODE.

# GROUPING TAXA BY GREPING FROM TXT FILE LIST
#	taxa_to_color <- read.table("/global/homes/r/rmbowers/databases_CogKeggPfam_master_tabs_GenomeTREELookuptabs_etc/genome_lookups_mostly_for_Trees_and_clusterMemberFAAs/taxa_to_color_full_phyla_colors_091316.txt", sep="\t", comment.char="")
	
	# if you want to color by whole clade
		# caveat - if clade is incorrectly named, i.e. 1 branch in wrong position, the mrca will be id'ed and colored, coloring too much
	#taxa_to_color <- as.character(taxa_to_color$V1)

	taxa_to_color <- read.table(args[3], header=FALSE, sep="\t")
	#taxa_to_color <- read.table("/projectb/scratch/rmbowers/DBASES_USEFUL_TO_NEW_LIFE_PROJECTS/arcbac27October2017_genomes_fna_faa/arcbac_taxa_color_names_rnaPolTree.txt", header=FALSE, sep="\t")
	#taxa_to_color <- read.table("/projectb/scratch/rmbowers/DBASES_USEFUL_TO_NEW_LIFE_PROJECTS/arcbac27October2017_genomes_fna_faa/arcbac_taxa_color_names_UNI56_Tree.txt", header=FALSE, sep="\t")
	taxa_to_color.tip <- lapply(taxa_to_color$V2, function(x) {tree_rerooted_rerooted$tip.label[grepl(x, tree_rerooted_rerooted$tip.label, ignore.case=FALSE)]})
	names(taxa_to_color.tip) <- paste(taxa_to_color$V1, taxa_to_color$V2, sep=";")
	taxa_to_color.mrca <- lapply(taxa_to_color.tip, function(x) {MRCA(tree_rerooted_rerooted, tip=x)})
	groups = list()
	groups <- lapply(taxa_to_color.mrca, function(x) {getDescendants(tree_rerooted_rerooted, node=x)})

	# if you want to color by grepping from names
		# caveat - only coloring those taxa with the name found in the grep() command
	#groups = list()
	#for(i in seq_along(taxa_to_color$V1)){groups[[i]] <- grep(taxa_to_color[[i, 1]], tree_rerooted_rerooted$tip.label)}
	#names(groups) <- taxa_to_color$V1
	#groups=groups[lapply(groups,length)>5]
	#taxa_to_color <- as.data.frame(taxa_to_color[match(names(groups), taxa_to_color$V1), ])
	#colnames(taxa_to_color) <- c("V1","V2")
	#taxa_to_color$V1 <- droplevels(as.factor(taxa_to_color$V1))

	tree_grouped <- groupOTU(tree_rerooted_rerooted, groups)
	#tree_grouped$tip.label <- sub(".*\\|", "", tree_grouped$tip.label)

# Coloring based on 'taxa_to_color' df, have not been able to figure out how to name colors/tree branches with correct order so they match the legend# https://www.r-bloggers.com/the-paul-tol-21-color-salute/
	tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
	colors <- tol21rainbow
	pal <- colorRampPalette(colors)
	branch.cols <- pal(nlevels(as.factor(names(groups))))

# place data into pseq object
	pseq_object <- phyloseq(otu_table(OTU, taxa_are_rows = TRUE), phy_tree(tree_grouped))

	sampleType_colors <- c("Samples_added_to_Tree" = "red")


if(args[4] == "linear") {


	file="tree_CPR_outgroup_linear.pdf"
	#pdf(file, width=200, height=200, pointsize=48) # circular
	pdf(file, width=120, height=500, pointsize=24) # linear
	tree_p <- ggtree(tree_grouped, ladderize=TRUE, size=5, alpha=1, aes(color=group)) +
		scale_color_manual(name="Phyla", values=c("gray70",  branch.cols)) +
		geom_tiplab(size=8, hjust= -0.2, show.legend=F) + #linear
		#geom_tiplab(size=8, aes(angle=angle), hjust= -0.2, show.legend=F) + #angled
		#geom_text(size=20, aes(label=node)) +
		geom_rootpoint() +
		#geom_text(size=20, aes(label=node)) +
		#ggtitle(z) +
		#theme(plot.title = element_text(size=20, face="bold.italic", hjust = 0.5)) +
		#scale_color_manual(values=sampleType_colors) +
		theme(legend.position="right") +
		theme(legend.title=element_blank(), legend.text=element_text(size=150, face="bold")) +
		guides(colour = guide_legend(override.aes = list(size=150), ncol=1)) +
		theme(panel.background = element_rect(fill="transparent", colour=NA),
					panel.grid.minor=element_blank(),
					panel.grid.major=element_blank(),
					plot.background=element_rect(fill="transparent", colour=NA)) +
		ggplot2::xlim(0, 5)

} else if(args[4] == "circular") {

	file="tree_CPR_outgroup_circular.pdf"
	pdf(file, width=400, height=400, pointsize=48) # circular
	#pdf(file, width=120, height=500, pointsize=24) # linear
	tree_p <- ggtree(tree_grouped, layout="circular", size=8, alpha=1, aes(color=group)) +
		scale_color_manual(name="Phyla", values=c("gray70",  branch.cols)) +
		#geom_tiplab(size=8, hjust= -0.2, show.legend=F) + #linear
		geom_tiplab(size=8, aes(angle=angle), hjust= -0.2, show.legend=F) + #angled
		#geom_text(size=20, aes(label=node)) +
		geom_rootpoint() +
		#geom_text(size=20, aes(label=node)) +
		#ggtitle(z) +
		#theme(plot.title = element_text(size=20, face="bold.italic", hjust = 0.5)) +
		#scale_color_manual(values=sampleType_colors) +
		theme(legend.position="right") +
		theme(legend.title=element_blank(), legend.text=element_text(size=150, face="bold")) +
		guides(colour = guide_legend(override.aes = list(size=150), ncol=1)) +
		theme(panel.background = element_rect(fill="transparent", colour=NA),
					panel.grid.minor=element_blank(),
					panel.grid.major=element_blank(),
					plot.background=element_rect(fill="transparent", colour=NA))}

# extracting data from pseq_object, so can plot in ggtree instead of phyloseq
	dd <- psmelt(pseq_object)
	# only want to plot points with abundances gtr than 0
	dd_gtr0 <- dd[dd$Abundance > 0, ]
	tree_p_temp_data <- tree_p$data
	# use sapply here b/c you want to perform the gsub function on a list, but then return a vector
	#tree_p_temp_data$label <- sapply(tree_p_temp_data$label, function(x) { gsub("\\'", "", x)})
	data_tree_p <- merge(tree_p_temp_data, dd_gtr0, by.x="label", by.y="OTU", all=TRUE)
	spacing <- 0.05
	idx <- with(data_tree_p, sapply(table(node)[unique(node)], function(i) 1:i)) %>% unlist
	hjust <- spacing * idx * max(data_tree_p$x)
	data_tree_p$xdodge <- data_tree_p$x + hjust
	data_tree_p$Sample <- as.factor(data_tree_p$Sample)
	data_tree_p <- na.omit(data_tree_p)
	

	### ADD ADDITIONAL METADATA TO DATA EXTRACTED FROM PHYLOSEQ if needed - not needed here
	# data_tree_p <- merge(data_tree_p, mapping_file[,c(2,5)], by.x="label", by.y="alternate_unique_id", all.x=TRUE)

## seems that you need the shape="" in order for fill="" to work correctly
tree_plot <- tree_p + geom_point(data=data_tree_p, aes(x=xdodge, size=Abundance, fill=Sample), alpha=1, shape=21) + theme(legend.position="right") +
	#geom_point(data=data_tree_p_16SNeg, aes(x=xdodge, size=Abundance, fill=Sample), alpha=0.2, shape=21) +
	#geom_point(data=data_tree_p_mtg, aes(x=xdodge, size=Abundance, fill=Sample), alpha=0.8, shape=21) + 
	scale_size(range=c(30)) +
	#scale_size_manual(values=c(1)) +
	scale_fill_manual(name="Sample Type", values=sampleType_colors) +
	#scale_shape_manual(name="Sample Type", values=c(23,21,23,21,23,21)) +
	#guides(shape=FALSE, fill=guide_legend(override.aes=list(shape=c(23,21,22))))
	guides(fill=guide_legend(override.aes = list(size=100), ncol=1)) 
	#guides(shape=guide_legend(override.aes = list(size=100))) 
	# unclustered point size
	# clustered point size
	#scale_size(range=c(50,120))	
# Cannot collapse root while also placing bubbles on tips. Tip placement will not correspond to bubbles 
#tree_p <- ggtree::collapse(tree_p, root_nodeID_rerootID) + geom_cladelabel(fontsize=40, offset=0.05, node=root_nodeID_rerootID, "Root")
tree_plot
dev.off()








