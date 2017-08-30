#' plotLmerMeans
#'
#' Given a model and a variable of iterest (and an interacting variable if required), code predicts new values based on reference values (first level of factors, any value of continuous variable specified by fun), and plots them.
#' @param mod.for.plot Model fit by lme4
#' @param do_plot If TRUE, plot the results. If FALSE, provide a table with the data for plotting.
#' @param seMultiplier Multiplier for plotting error bars. Defaults to 1.96 which will provide 95% confidence intervals.
#' @param with_ranefs Defaults to FALSE. If TRUE, will include some error in the intercept associated with the random effects. I don't think you'll ever really want to do this, so stick with FALSE.
#' @param var.of.interest Variable you want to plot
#' @param var.levels If you want to plot only some levels of var.of.interest, add them here
#' @param var.interaction Variable interacting with var.of.interest that you wish to plot (leave blank if you want main effects)
#' @param var.interaction.levels If you want to plot only some levels of var.interaction, add them here
#' @param fun The function that determines the baseline level of continuous variables on which to base predictions (e.g., max, min, median)
#' @param length If var.of.interest or var.interaction is a continuous variable, give the length of the sequence you want to predict over (between 50 and 200 usually suffice)
#' @param ylabel,xlabel,mainlabel The x, y and main labels for figure
#' @param meancol,CIcol Colours for main points/lines and error bars
#' @param point.size,point_type,line.width Point size (cex), Point type (pch), Line width (lwd)
#' @param transform_fun If you wish to back-transform data, provide function here e.g. exp
#' @param with_CI Whether to include error bars
#' @param CI_fill For continuous variables, if TRUE, fill the confidence intervals with colour. If FALSE, add confidence intervals as lines.
#' @param return_quietly If FALSE, provide table of data for plotting
#' @param with_labs Whether to include labels for categorical variables on the x axis
#' @param with_ticks Whether to include tick marks on y axis
#' @param with_xaxis Whether to include an x axis at all
#' @param lab_cex character expansion for labels
#' @param axis_cex character expansion for axis
#' @param ylims If you want to specify the upper and lower limits on the yaxis, do so here. Otherwise limits will be automatically calculated.
#' @examples
#' plotLmerMeans(m1, var.of.interest = "LandUse", var.interaction = "BeforeAfter")
#' @note This is not a pretty function and has grown over the years, so let me know if there are any issues. The response variable must be identifiable within the dataframe (i.e., log(data) will not work - you must have another variable in the datframe where logged data exists)
#' @author Adriana De Palma
#' @export
#' @import stats

plotLmerMeansRepeat <- function(mod.for.plot, do_plot = TRUE, seMultiplier = 1.96, with_ranefs = FALSE,
                        var.of.interest = NULL, var.levels = NULL, var.interaction = NULL, var.interaction.levels = NULL,
                        fun = min,length = 200,ylabel = "", xlabel = "", mainlabel = "",
                        CIcol = "light blue", meancol ="#1F78B4",
                        point.size = 2,point_type = 16,line.width = 5,
                        transform_fun,
                        with_CI = TRUE, CI_fill = FALSE,return_quietly = TRUE,
                        with_labs = TRUE,with_ticks = TRUE, with_xaxis = TRUE,lab_cex = 1, axis_cex = 1,ylims){
  
  
  require(gstat, quietly = TRUE)
  require(lattice, quietly = TRUE)
  require(mgcv, quietly = TRUE)
  require(scatterplot3d, quietly = TRUE)
  
  
  random.vars<-mod.for.plot@cnms ## extract the random variables from the model
  
  response.variable<-names(mod.for.plot@frame)[1] ## extract the name of the response variable
  other.variables<-names(mod.for.plot@frame)[2:(length(names(mod.for.plot@frame))-length(random.vars))]  ## extract the name of the other variables (the random effects appear at the end of this list (not including residual) hence -length(random.actual))
  
  ## if there are other terms in the model
  term.levs<-list()
  for(i in 1:length(other.variables)){ ## for all terms in the model
    
    if(is.factor(mod.for.plot@frame[,other.variables[i]])==TRUE){ ## if it is a factor
      term.levs[[i]]<-levels(mod.for.plot@frame[,other.variables[i]]) ## get its levels from the model frame
    }
    else{ ## if it is a continuous variable
      if(grepl(var.of.interest,other.variables[i])){ ## if the term is the variable of interest
        term.levs[[i]]<-seq(min(mod.for.plot@frame[,other.variables[i]]),max(mod.for.plot@frame[,other.variables[i]]),length = length)
        # then get the full sequence
      }else{
        if(!is.null(var.interaction)){
          if(grepl(var.interaction,other.variables[i])){
            term.levs[[i]]<-seq(min(mod.for.plot@frame[,other.variables[i]]),max(mod.for.plot@frame[,other.variables[i]]),length = length)
          }else{
            term.levs[[i]]<-fun(mod.for.plot@frame[,other.variables[i]])
          }
        }else{
          term.levs[[i]]<-fun(mod.for.plot@frame[,other.variables[i]])
        } ## if it is not a variable of interest, just take the fun of the sequence
        
      }
    }
    
  } # close the loop
  
  newdat<-expand.grid(term.levs)
  names(newdat)<-other.variables
  newdat[,response.variable]<-rep(0,nrow(newdat))
  
  ## minimise the predictions to what we're after
  newdat.default<-newdat
  
  if(!is.null(var.levels)){
    newdat.default<-newdat[newdat[,var.of.interest]%in%var.levels,] ## only use levels specified
  }
  
  vterms<-other.variables[other.variables!=var.of.interest] ## remove the factor of interest
  
  if(!is.null(var.interaction)){ ## if there is an interacting term to plot
    terms.def<-vterms[vterms!=var.interaction]
    if(!is.null(var.interaction.levels)){ ## if we want all the terms for this interaction
      newdat.default<-newdat.default[newdat.default[,var.interaction]%in%var.interaction.levels,]
    }
  }else{
    terms.def<-vterms
  }
  if(length(terms.def)!=0){## if there is more than just one explanatory variable in the model
    for(i in 1:length(terms.def)){ ## for each additional explanatory variable
      if(is.factor(mod.for.plot@frame[,terms.def[i]])==TRUE){ ## get the default values of the remaining terms
        iterm<-terms.def[i]
        default.intercept<-levels(mod.for.plot@frame[,iterm])[length(levels(mod.for.plot@frame[,iterm]))]
        newdat.default<-newdat.default[newdat.default[,iterm]==default.intercept,]
      }
    }
  } ## subsets the predictions so it's only using the values for the reference levels of the other terms. Else leave newdat as is
  
  ## minimises the mm to make sure that there's no rank deficiency
  coef.names<-names(fixef(mod.for.plot))
  mm<-model.matrix(terms(mod.for.plot), newdat.default)
  original.mm.names<-dimnames(mm)[[2]]
  if(!missing(var.interaction)){
    if(length(coef.names)!=length(original.mm.names)){
      mm<-mm[,dimnames(mm)[[2]]%in%coef.names]
    }
  }
  
  
  newdat.default[,response.variable]<-mm %*% fixef(mod.for.plot) ## point predictions
  
  pvar1 <- diag(mm %*% tcrossprod(as.matrix(vcov(mod.for.plot)),mm))## variance without random effects
  ran.sum<-0
  for(i in 1:length(random.vars)){
    ran.sum<-ran.sum+VarCorr(mod.for.plot)[[i]][1]
  }
  tvar1<-pvar1+ran.sum ## adding up the fixed effects variance and random variances - do this for all random effects
  newdat.default<-data.frame(
    newdat.default
    , plo = newdat.default[,response.variable]-sqrt(pvar1) * seMultiplier
    , phi = newdat.default[,response.variable]+sqrt(pvar1) * seMultiplier
    , tlo = newdat.default[,response.variable]-sqrt(tvar1) * seMultiplier
    , thi = newdat.default[,response.variable]+sqrt(tvar1) * seMultiplier
  )
  if(missing(transform_fun)==FALSE){
    newdat.default[,response.variable]<-transform_fun(newdat.default[,response.variable])
    newdat.default$plo<-transform_fun(newdat.default$plo)
    newdat.default$phi<-transform_fun(newdat.default$phi)
    newdat.default$tlo<-transform_fun(newdat.default$tlo)
    newdat.default$thi<-transform_fun(newdat.default$thi)
  }
  
  
  
  #minimise newdat.default so that we don't have spurious values for coefficients that weren't estimated
  
  if(!missing(var.interaction)){
    if(is.factor(mod.for.plot@frame[,var.interaction])&is.factor(mod.for.plot@frame[,var.of.interest])){
      if(length(coef.names)!=length(original.mm.names)){
        
        t<-as.data.frame(table(mod.for.plot@frame[,var.of.interest],mod.for.plot@frame[,var.interaction]))
        names(t)<-c(var.of.interest,var.interaction,"Freq")
        
        alllevs<-paste(newdat.default[,var.of.interest],newdat.default[,var.interaction])
        t$levs<-paste(t[,var.of.interest],t[,var.interaction])
        
        t2<-t[match(alllevs,t$levs),]
        
        stopifnot(t2$levs%in%alllevs)
        
        newdat.default<-newdat.default[t2$Freq!=0,]
      }
    }
  }
  
  
  if (do_plot == TRUE){ ## otherwise plot the results for fac of interest
    if(!missing(meancol)){
      if(missing(CIcol)){
        CIcol<-meancol
      }
    }
    if(!missing(CIcol)){
      if(missing(meancol)){
        meancol<-CIcol
      }
    }
    
    ## set up the reusable plotting functions
    
    points.plot<-function(preds,levs = "lev", ranef = TRUE,labs = TRUE){
      if(length(meancol)==1){
        meancol<-rep(meancol,nrow(preds))
      }
      else{
        if(length(meancol)<nrow(preds)){
          print("Warning: not enough plotting colours specified - using first colour only")
          meancol<-rep(meancol[1],nrow(preds))
        }
        if(length(meancol)>nrow(preds)){
          print("Warning: too many plotting colours specified - using in order given ")
        }
      }
      if(length(CIcol)==1){
        CIcol<-rep(CIcol,nrow(preds))
      }
      else{
        if(length(CIcol)<nrow(preds)){
          print("Warning: not enough plotting colours specified - using first colour only")
          CIcol<-rep(CIcol[1],nrow(preds))
        }
        if(length(CIcol)>nrow(preds)){
          print("Warning: too many plotting colours specified - using in order given ")
        }
      }
      
      if(length(point_type)==1){
        point_type<-rep(point_type,nrow(preds))
      }
      else{
        if(length(point_type)<nrow(preds)){
          print("Warning: not enough plotting characters specified - using first character only")
          point_type<-rep(point_type[1],nrow(preds))
        }
        if(length(point_type)>nrow(preds)){
          print("Warning: too many plotting characters specified - using in order given ")
        }
      }
      if(labs == TRUE){
        for(i in 1:length(preds[,levs])){
          mtext(preds[i,levs],side = 1, line = linenum, at = i, adj = 0.5)
        }
      }
      if(ranef == TRUE){
        for(i in 1:nrow(preds)){
          arrows(i,preds$tlo[i],i,preds$thi[i],length = 0,col = CIcol[i],lwd = line.width)
          points(i,preds[i,response.variable],cex = point.size, col = meancol[i],pch = point_type[i])
        }
      }
      if(ranef == FALSE){
        for(i in 1:nrow(preds)){
          arrows(i,preds$plo[i],i,preds$phi[i],length = 0,col = CIcol[i],lwd = line.width)
          points(i,preds[i,response.variable],cex = point.size, col = meancol[i],pch = point_type[i])
        }
      }
    }
    
    lines.plot<-function(preds,vars,ranef = FALSE, with_CI = TRUE){
      if(length(vars)>2){
        print("too many variables specified")
      }
      else{
        if(length(vars)==0){
          print("no variables specified")
        }
      }
      for(i in 1:length(vars)){
        if(is.factor(preds[,vars[i]])==TRUE){
          fac.var<-vars[i]
          n.lines<-length(levels(preds[,vars[i]]))
          fac.levs<-levels(preds[,vars[i]])
        }
        else{
          cont.var<-vars[i]
        }
      }
      if(length(vars)==1){ ## if there's no interaction
        if(length(meancol)>1){
          print("Warning: too many plotting colours specified - only the first used ")
          meancol<-meancol[1]
        }
        if(ranef==TRUE){
          if(with_CI == TRUE){
            if(CI_fill == TRUE){
              polygon(c(preds[,cont.var],rev(preds[,cont.var])), c(preds$thi,rev(preds$tlo)),col = CIcol[1],border = CIcol[1],lwd = 2)
            }
            else{
              lines(preds[,cont.var],preds$thi,col = CIcol[1],lwd = 2, lty = 2)
              lines(preds[,cont.var],preds$tlo,col = CIcol[1],lwd = 2, lty = 2)
            }
          }
          lines(preds[,cont.var],preds[,response.variable],col = meancol[1],lwd = line.width)
        }
        else{
          if(with_CI == TRUE){
            if(CI_fill == TRUE){
              polygon(c(preds[,cont.var],rev(preds[,cont.var])), c(preds$phi,rev(preds$plo)),col = CIcol[1],border = CIcol[1],lwd = 2)
            }
            else{
              lines(preds[,cont.var],preds$phi,col = CIcol[1],lwd = 2, lty = 2)
              lines(preds[,cont.var],preds$plo,col = CIcol[1],lwd = 2, lty = 2)            }
          }
          lines(preds[,cont.var],preds[,response.variable],col = meancol[1],lwd = line.width)
        }
      }
      else{ ## otherise if there is an interaction
        if(length(meancol)==1){
          meancol<-rep(meancol,n.lines)
        }
        else{
          if(length(meancol)<n.lines){
            print("Warning: not enough plotting colours specified - using first colour only")
            meancol<-rep(meancol[1],n.lines)
          }
          if(length(meancol)>n.lines){
            print("Warning: too many plotting colours specified - using in order given ")
          }
        }
        if(length(CIcol)==1){
          CIcol<-rep(CIcol,n.lines)
        }
        else{
          if(length(CIcol)<n.lines){
            print("Warning: not enough plotting colours specified - using first colour only")
            CIcol<-rep(CIcol[1],n.lines)
          }
          if(length(CIcol)>n.lines){
            print("Warning: too many plotting colours specified - using in order given")
          }
        }
        
        for(i in 1:n.lines){
          preds.i<-preds[preds[,fac.var]==fac.levs[i],] #subset for each level of factor
          if(ranef==TRUE){
            if(with_CI==TRUE){
              if(CI_fill == TRUE){
                polygon(c(preds.i[,cont.var],rev(preds.i[,cont.var])), c(preds.i$thi,rev(preds.i$tlo)),col = CIcol[i],border = CIcol[i],lwd = 2,lty = 1)
              }
              else{
                lines(preds.i[,cont.var],preds.i$thi,col = CIcol[i],lwd = 2, lty = 2)
                lines(preds.i[,cont.var],preds.i$tlo,col = CIcol[i],lwd = 2, lty = 2)
              }
            }
            lines(preds.i[,cont.var],preds.i[,response.variable],col = meancol[i],lwd = line.width)
          }
          else{
            if(with_CI==TRUE){
              if(CI_fill == TRUE){
                polygon(c(preds.i[,cont.var],rev(preds.i[,cont.var])), c(preds.i$phi,rev(preds.i$plo)),col = CIcol[i],border = CIcol[i],lwd = 2,lty = 1)
              }
              else{
                lines(preds.i[,cont.var],preds.i$phi,col = CIcol[i],lwd = 2, lty = 2)
                lines(preds.i[,cont.var],preds.i$plo,col = CIcol[i],lwd = 2, lty = 2)
              }
            }
            lines(preds.i[,cont.var],preds.i[,response.variable],col = meancol[i],lwd = line.width)
          }
        }
      }
    }
    
    
    
    
    
    ##make sure that var.of.interest and var.interaction are factors again if they need to be
    if(is.factor(mod.for.plot@frame[,var.of.interest])==TRUE){
      newdat.default[,var.of.interest]<-factor(newdat.default[,var.of.interest])
      newdat.default<-newdat.default[order(newdat.default[,var.of.interest]),]
    }
    else{
      newdat.default<-newdat.default[order(newdat.default[,var.of.interest]),]
    }
    if(!missing(var.interaction)){
      if(is.factor(mod.for.plot@frame[,var.interaction])==TRUE){
        newdat.default[,var.interaction]<-factor(newdat.default[,var.interaction])
      }
      else{
        newdat.default<-newdat.default[order(newdat.default[,var.interaction]),]
      }
    }
    
    
    
    ## set up plot
    
    if(!missing(ylims)){
      ylimlo <- min(ylims)
      ylimup<-max(ylims)
    }else{
      if(with_ranefs == TRUE){
        ylimlo<-min(newdat.default$tlo)-diff(range(newdat.default$tlo))
        ylimup<-max(newdat.default$thi)+diff(range(newdat.default$thi))
      }
      else{
        ylimlo<-min(newdat.default$plo)-diff(range(newdat.default$plo))
        ylimup<-max(newdat.default$phi)+diff(range(newdat.default$phi))
      }
    }
    
    
    
    
    if(is.factor(newdat.default[,var.of.interest])==TRUE){ ## if var of interest is a factor
      if(missing(var.interaction)){ ## and there are no interactions
        xlimlo<-0
        xlimup<-nrow(newdat.default)+1
        plot(NULL,ylim = c(ylimlo,ylimup), xlim = c(xlimlo,xlimup), ylab = ylabel,xlab = xlabel,axes = FALSE,main = mainlabel, cex.axis =axis_cex, cex.lab =lab_cex)
        axis(2,labels = TRUE)
        if(with_xaxis == TRUE){
          if(with_ticks == FALSE){
            axis(1,labels = FALSE,at = c(1:(xlimup-1)),lwd.ticks = -1)
          } else{
            axis(1,labels = FALSE,at = c(1:(xlimup-1)))
          }
        }
        
        newdat.default$lev<-paste(newdat.default[,var.of.interest])
        if(max(nchar(newdat.default$lev))>6){
          linenum<-5
        }
        else{
          linenum<-2
        }
        points.plot(preds = newdat.default,levs = "lev",ranef = with_ranefs,labs = with_labs)
      }
      else{ ## if var.of.interest is a factor and there is an interacting variable
        if(is.factor(newdat.default[,var.interaction])==TRUE){     ## if the interacting variable is a factor
          newdat.default$lev<-paste(newdat.default[,var.of.interest],newdat.default[,var.interaction],sep = ".")
          if(max(nchar(newdat.default$lev))>6){
            linenum<-5
          }
          else{
            linenum<-2
          }
          xlimlo<-0
          xlimup<-nrow(newdat.default)+1
          plot(NULL,ylim = c(ylimlo,ylimup), xlim = c(xlimlo,xlimup), ylab = ylabel,xlab = xlabel,axes = FALSE,main = mainlabel, cex.lab =lab_cex, cex.axis =axis_cex)
          axis(2,labels = TRUE)
          if(with_xaxis == TRUE){
            if(with_ticks == FALSE){
              axis(1,labels = FALSE,at = c(1:(xlimup-1)),lwd.ticks = -1)
            } else{
              axis(1,labels = FALSE,at = c(1:(xlimup-1)))
            }
          }
          points.plot(preds =newdat.default,levs = "lev",ranef = with_ranefs,labs = with_labs)
        }
        else{ # if the interacting variable is continuous
          xlimup<-max(newdat.default[,var.interaction])
          xlimlo<-min(newdat.default[,var.interaction])
          plot(NULL,ylim = c(ylimlo,ylimup), xlim = c(xlimlo,xlimup), ylab = ylabel,xlab = xlabel,axes = TRUE,main = mainlabel, cex.axis =axis_cex, cex.lab =lab_cex)
          lines.plot(newdat.default,vars = c(var.interaction,var.of.interest),ranef = with_ranefs,with_CI = with_CI)
        }
        
      }
    }
    else{ ## if the var.of.interest is continuous
      if(missing(var.interaction)){ ## if theres' no interaction
        xlimup<-max(newdat.default[,var.of.interest])
        xlimlo<-min(newdat.default[,var.of.interest])
        plot(NULL,ylim = c(ylimlo,ylimup), xlim = c(xlimlo,xlimup), ylab = ylabel,xlab = xlabel,axes = TRUE,main = mainlabel, cex.lab =lab_cex, cex.axis =axis_cex)
        lines.plot(preds = newdat.default,var = var.of.interest,ranef = with_ranefs,with_CI = with_CI)
      }
      else{ ## if there is an interacting variable
        if(is.factor(newdat.default[,var.interaction])){ ## and it is a factor
          xlimup<-max(newdat.default[,var.of.interest])
          xlimlo<-min(newdat.default[,var.of.interest])
          plot(NULL,ylim = c(ylimlo,ylimup), xlim = c(xlimlo,xlimup), ylab = ylabel,xlab = xlabel,axes = TRUE,main = mainlabel, cex.axis =axis_cex, cex.lab =lab_cex )
          lines.plot(preds = newdat.default,vars = c(var.of.interest,var.interaction),ranef = with_ranefs, with_CI = with_CI)
        }
        else{
          require(scatterplot3d,quietly = TRUE)
          scatterplot3d(newdat.default[,var.of.interest],newdat.default[,var.interaction],newdat.default[,response.variable],xlab = var.of.interest,zlab = response.variable,ylab = var.interaction,main = mainlabel)
        }
      }
    }
    if(return_quietly == FALSE){
      return(newdat.default)
    }
  } else{
    return(newdat.default) ## return dataframe of predictions
  }
}
