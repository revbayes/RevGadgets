################################################################################
#
# @brief Plotting the output of a episodic diversification rate analysis with mass-extinction events.
#
# @date Last modified: 2016-09-05
# @author Sebastian Hoehna
# @version 1.0
# @since 2018-12-06, version 1.0.0
#
# @param    output            list          The processed output for plotting.
# @param    fig.types         character     Which aspects of the model to visualize. See details for a complete description.
# @param    xlab              character     The label of the x-axis. By default, millions of years.
# @param    col               character     Colors used for printing. Must be of same length as fig.types.
# @param    col.alpha         numeric       Alpha channel parameter for credible intervals.
# @param    xaxt              character     The type of x-axis to plot. By default, no x-axis is plotted (recommended).
# @param    yaxt              character     The type of y-axis to plot.
# @param    pch               integer       The type of points to draw (if points are drawn).
# @param    ...                             Parameters delegated to various plotting functions.
#
#
################################################################################

rev.plot.pop.size = function(output,
                             fig.types=c("population size", "population size shift prob"),
                             time=100,
                             xlab="years ago",col=NULL,col.alpha=50,
                             xaxt="n",yaxt="s",pch=19,
                             use.smoothing=FALSE,
                             ...){

#fig.types=c("population size", "population size shift prob")

#    # Check that fig type is valid
#    validFigTypes <- c("speciation rate","speciation shift times","speciation Bayes factors",
#                       "extinction rate","extinction shift times","extinction Bayes factors",
#                       "fossilization rate","fossilization shift times","fossilization Bayes factors",
#                       "net-diversification rate","relative-extinction rate",
#                       "mass extinction times","mass extinction Bayes factors")
#    invalidFigTypes <- fig.types[!fig.types %in% validFigTypes]
#
#    if ( length( invalidFigTypes ) > 0 ) {
#        stop("\nThe following figure types are invalid: ",paste(invalidFigTypes,collapse=", "),".",
#             "\nValid options are: ",paste(validFigTypes,collapse=", "),".")
#    }

    # Make color vector
    if ( is.null(col) ) {
        col <- c("population size"="#984EA3",
                 "population size shift times"="#984EA3",
                 "population size shift prob"="#984EA3",
                 "population size Bayes factors"="#984EA3")
    } else {
        names(col) <- fig.types
    }

    # Compute the axes
    tree_age <- time
#    tree_age <- max(node.depth.edgelength(tree))
    num_intervals <- length(output$intervals)-1
    plot_at <- 0:num_intervals
    interval_size <- tree_age/num_intervals
    labels <- pretty(c(0,tree_age))
    labels_at <- num_intervals - (labels / interval_size)
    ages <- seq(0,tree_age,length.out=num_intervals+1)

    for( type in fig.types ) {

        if ( grepl("times",type) ) {

            thisOutput <- output[[type]]
            meanThisOutput <- colMeans(thisOutput)
            criticalPP <- output[[grep(strsplit(type," ")[[1]][1],grep("CriticalPosteriorProbabilities",names(output),value=TRUE),value=TRUE)]]

            if(plot.tree){
                plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,treeAge))
                par(new=TRUE)
            }

            barplot(meanThisOutput,space=0,xaxt=xaxt,col=col[type],border=col[type],main=type,ylab="posterior probability",xlab=xlab,ylim=c(0,1),...)
            abline(h=criticalPP,lty=2,...)
            axis(4,at=criticalPP,labels=2*log(output$criticalBayesFactors),las=1,tick=FALSE,line=-0.5)
            axis(1,at=labelsAt,labels=labels)
            box()

        } else if ( grepl("Bayes factors",type) ) {

            thisOutput <- output[[type]]
            ylim <- range(c(thisOutput,-10,10),finite=TRUE)

            if(plot.tree){
                plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,tree_age))
                par(new=TRUE)
            }
            plot(x=plotAt[-1]-diff(plotAt[1:2])/2,y=thisOutput,type="p",xaxt=xaxt,col=col[type],ylab="Bayes factors",main=type,xlab=xlab,ylim=ylim,xlim=range(plotAt),pch=pch,...)
            abline(h=2 * log(output$criticalBayesFactors),lty=2,...)
            axis(4,at=2 * log(output$criticalBayesFactors),las=1,tick=FALSE,line=-0.5)
            axis(1,at=labelsAt,labels=labels)

        } else if ( grepl("shift prob",type) ) {

            this_output <- output[[type]]
            ylim <- range(c(this_output,0,1),finite=TRUE)

#            if(plot.tree){
#                plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,tree_age))
#                par(new=TRUE)
#            }
            plot(x=plot_at[-1]-diff(plot_at[1:2])/2,y=this_output,type="p",xaxt=xaxt,col=col[type],ylab="Posterior Probability",main=type,xlab=xlab,ylim=ylim,xlim=range(plot_at),pch=pch,cex=0.1,...)
#            abline(h=2 * log(output$criticalBayesFactors),lty=2,...)
#            axis(4,at=2 * log(output$criticalBayesFactors),las=1,tick=FALSE,line=-0.5)
            axis(1,at=labels_at,labels=labels)

        } else {

            this_output <- output[[type]]
            mean_this_output <- colMeans(this_output)      
            quantiles_this_output <- apply(this_output,2,quantile,prob=c(0.025,0.975))
            
            if ( use.smoothing == TRUE ) {
                y_unsmoothed <- c(mean_this_output[1],mean_this_output)
                y_smoothed <- (lowess(x=1:length(y_unsmoothed),y_unsmoothed, f=0.2))$y
                    
                y_unsmoothed_1 <- c(quantiles_this_output[1,1],quantiles_this_output[1,])
                y_unsmoothed_2 <- rev(c(quantiles_this_output[2,1],quantiles_this_output[2,]))
                lo_1 <- lowess(x=1:length(y_unsmoothed_1),y_unsmoothed_1, f=0.2)
                y_smoothed_1 <- lo_1$y
                lo_2 <- lowess(x=1:length(y_unsmoothed_2),y_unsmoothed_2, f=0.2)
                y_smoothed_2 <- lo_2$y
                y_smoothed_ci <- c(y_smoothed_1,y_smoothed_2)
            } else {
                y_smoothed <- c(mean_this_output[1],mean_this_output)
                y_smoothed_ci <- c(c(quantiles_this_output[1,1],quantiles_this_output[1,]),rev(c(quantiles_this_output[2,1],quantiles_this_output[2,])))
            }

            
            # find the limits of the y-axis
            ylim <- c(min(0,quantiles_this_output),max(quantiles_this_output))

#            if (plot.tree) {
#                plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,tree_age))
#                par(new=TRUE)
#            }

            
#            plot(x=plot_at,y=c(mean_this_output[1],mean_this_output),type="l",ylim=ylim,xaxt=xaxt,col=col[type],ylab="Population size",main=type,xlab=xlab,...)
#            polygon(x=c(0:ncol(quantiles_this_output),ncol(quantiles_this_output):0),y=c(c(quantiles_this_output[1,1],quantiles_this_output[1,]),rev(c(quantiles_this_output[2,1],quantiles_this_output[2,]))),border=NA,col=paste(col[type],col.alpha,sep=""))
            
            
            plot(x=plot_at,y=y_smoothed,type="l",ylim=ylim,xaxt=xaxt,col=col[type],ylab="Population size",main=type,xlab=xlab,...)
            polygon(x=c(0:ncol(quantiles_this_output),ncol(quantiles_this_output):0),y=y_smoothed_ci,border=NA,col=paste(col[type],col.alpha,sep=""))

            axis(1,at=labels_at,labels=labels)
            
      
        }

    }

}



################################################################################
#
# @brief Processing the output of a "skyline" coalescent analysis.
#
# @date Last modified: 2018-12-06
# @author Sebastian Hoehna
# @version 1.0
# @since 2018-12-06, version 1.0.0
#
# @param    numExpectedRateChanges      numeric        The number of expected diversification-rate changes.
# @param    numExpectedMassExtinctions  numeric        The number of expected mass-extinction events.
# @param    burnin                      numeric        The fraction of the posterior samples to be discarded as burnin.
# @param    numIntervals                numeric        The number of discrete intervals in which to break the tree.
#
#
################################################################################

rev.process.pop.size = function(population_size_change_times_file="",population_size_file="", time=0, burnin=0.25, num_intervals=100){


    # Get the time of the tree and divide it into intervals
    intervals <- seq(0,time,length.out=num_intervals+1)

    process_pop_sizes <- rev.read.mcmc.output.rates.through.time(population_size_change_times_file, population_size_file, intervals, burnin)
    
    process_pop_size_shift_prob <- c()
    for ( i in 2:num_intervals ) {
        process_pop_size_shift_prob[i] <- mean(process_pop_sizes[,i] > process_pop_sizes[,i-1]) + 0.5 * mean(process_pop_sizes[,i] == process_pop_sizes[,i-1])
    }

    res <- list("population size" = process_pop_sizes,
                "population size shift prob" = process_pop_size_shift_prob,
                "intervals" = rev(intervals) )

    return(res)
}

