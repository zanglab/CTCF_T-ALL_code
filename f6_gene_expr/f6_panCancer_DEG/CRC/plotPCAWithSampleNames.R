# Modified DESeq2 plotPCA function with sample names and proportion of variance added.
# Sample names will be shown underneath each dot.
# The axis will display proportion of variance for each principal component.
# Tested using DESeq2 1.2.8, 1.6.2, and 1.8.1.
# The native DESeq2 plotPCA function switched from lattice to ggplot2 in version 1.5.11.

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=500,label_col=1)
{
	library(RColorBrewer)
	library(genefilter)
	library(lattice)
    #group_name_table = c("NCE2.0","NCE2.5a","NCE3.0","MEKi-NCE2.0","C-NCE2.5a","E-NCE2.5a","NCE2.5b")
	# pca
	rv = rowVars(assay(x))
	select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
	pca = prcomp(t(assay(x)[select,]))

	# proportion of variance
	variance = pca$sdev^2 / sum(pca$sdev^2)
	variance = round(variance, 3) * 100
	
	# sample names
	names = colnames(x)

	# factor of groups
	fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))

	# colors
	if( nlevels(fac) >= 10 )
		colors = rainbow(nlevels(fac))
	else if( nlevels(fac) >= 3 )
		colors = brewer.pal(nlevels(fac), "Set1")
	else
		colors = c(  "firebrick3","dodgerblue3" )
    
    #colors_set = brewer.pal(9, "Set1")
    #c(colors_set[1:3],colors_set[7],colors_set[5:6])
    #colors = colors_set[sort(x$condition[!duplicated(x$condition)])]
	# plot
	plt<- xyplot(
		PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.3,
		aspect = 1,layout=c(1, 1),
		#type='a'
		col = colors,
		#auto.key =
        # list(space = "right", points = FALSE, lines = TRUE, cex=0.25)
		xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.9),
		ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.9),
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...);
			#ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
		},
		#main = draw.key(
			key = list(
			    space='top',#x=0.92,y=0.3,
			    corner = c(1,0),  
				columns=2,text.width=0,pch=c(19,19),
				rect = list(col = colors),
				#text = list(group_name_table[sort(x$condition[!duplicated(x$condition)])]),
				text = list(levels(fac)),
				rep = FALSE,cex=0.8,between=0.3,size=2
				#,pch=10,padding.text=0.2,size=2
			#)
		)
		,lattice.options = list(legend.bbox="full") 
	)

    update(plt, par.settings = list(fontsize = list(text = 12, points = 8)))

}
