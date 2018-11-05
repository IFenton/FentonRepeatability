# code for plotting confusion matrices

conf_mat <- function(data.long, col1, col2 = col1, axis.col = NULL, axis1 = NULL, axis2 = NULL, key = TRUE, num1 = "full", num2 = "sample", spec.abb = NULL, abb.end = NULL, axes.same = TRUE, xlab = paste(col1), ylab = paste(col2), sp.list = "present", sp.exc = NULL, subset.col = NULL, subset.lev = NULL, palette = "matlab.like", grid = FALSE, grid.sp = 5, grid.col = "grey40") {
  # dataset in long form
  # col 1 / 2 - columns from the dataset 
  # axis 1 / 2 - if the columns are the same, what subset of the data are being used?
  # num 1 / 2 - should numbers be "full" for all samples, or just for that "sample"
  # spec.abb - the link between full species names and abbreviations, colnames "Species" "Abbreviations". If NULL, it assumes that the names in the table are the full names
  # abb.end - any abbreviations that should be moved to the end of the sort
  # axes.same - should the x and y axes have the same species on them? x will always have the full list, y can be shorter
  # sp.list - should only the species "present" in the subset of the data be used, or should the "full" list be used
  # sp.exc - species to be excluded from the analysis
  # subset.col / subset.lev - the column and level by which to subset the data
  require(caret) # for the confusion matrix
  require(colorRamps) # colours
  require(stringr)
  par.def <- par("fig", "mar", "mai")
  on.exit(par(par.def))
  if (key) {
    par(fig = c(0, 0.9, 0, 1))
  }
  par(mar = c(16, 16, 2.5, 2))
  # generate a species abbreviations table if it doesn't already exist
  if (is.null(spec.abb)) {
    spec.abb <- data.frame(Species = sort(as.character(unique(data.long[, col1]))), Abbreviations = sort(as.character(unique(data.long[, col1]))), stringsAsFactors = FALSE)
  }
  # create a dataframe to work on
  # if there are some species to exclude, e.g. "lost"
  if (!is.null(sp.exc)) {
    tmp.exc <- c(which(data.long[, col1] == sp.exc), which(data.long[, col2] == sp.exc))
    data <- data.long[-tmp.exc, ]
  } else {
    data <- data.long
  }
  # if we should only use a subset of the data
  if (!is.null(subset.col)) {
    if (is.null(subset.lev))
      stop("subset level is needed when subset.col is used")
    data <- data[data[,subset.col] == subset.lev, ]
  }
  
  
  # create the smaller dataframe for the analysis
  if (is.null(axis.col)) {
    data <- data.frame(ax1 = as.character(data[, col1]), ax2 = as.character(data[, col2]), stringsAsFactors = FALSE)
  } else {
    data <- data.frame(ax1 = as.character(data[data[,names(data) == axis.col] == axis1, names(data) == col1]), ax2 = as.character(data[data[,names(data) == axis.col] == axis2, names(data) == col2]), stringsAsFactors = FALSE)
  }
  
  # which subset of the species are needed for this analysis
  if (sp.list == "present") {
    spec.abb <- spec.abb[spec.abb$Abbreviation %in% unique(c(data$ax1, data$ax2)),]
    # reorder
    spec.abb <- spec.abb[c(which(!(spec.abb$Abbreviation %in% abb.end)), which(spec.abb$Abbreviation %in% abb.end)), ]
  } else if (sp.list == "full") {
    spec.abb <- spec.abb[c(which(!(spec.abb$Abbreviation %in% abb.end)), which(spec.abb$Abbreviation %in% abb.end)), ]
  }
  
    # species for axis 1 / axis 2
  if (axes.same) {
    spec.ax1 <- spec.ax2 <- spec.abb
  } else {
    # reorder so that the species in the shorter axis come first
    spec.abb <- spec.abb[c(which(spec.abb$Abbreviation %in% unique(data[, 2])), which(!(spec.abb$Abbreviation %in% unique(data[, 2])))), ]
    spec.ax1 <- spec.abb
    spec.ax2 <- spec.abb[spec.abb$Abbreviation %in% unique(data[, 2]),]
  }
  
  # the confusion matrix
  conf <- confusionMatrix(factor(data$ax2, levels = spec.abb$Abbreviation), factor(data$ax1, levels = spec.abb$Abbreviation))
  
  dim1 <- nrow(spec.ax1)
  dim2 <- nrow(spec.ax2)
  xlim1 <- -(1+1/(dim1 - 1))/dim1/2
  xlim2 <- (dim1*((1+1/(dim1 - 1)))/dim1) - (1+1/(dim1 - 1))/dim1/2
  ylim1 <- (dim2*((1+1/(dim1 - 1)))/dim1) - (1+1/(dim1 - 1))/dim1/2
  ylim2 <- -(1+1/(dim1 - 1))/dim1/2
  
  plot(1, type = "n", xlim = c(xlim1, xlim2), ylim = c(ylim1, ylim2), xlab = "", ylab = "", axes = FALSE, bty = "o", xaxs = "i", yaxs = "i")
  rect(xlim1, ylim2, xlim2, ylim1, col = "grey70")
  image(t(conf$table / rowSums(conf$table)), col =c("grey70", do.call(palette, list(100000))), axes = FALSE, add = TRUE)
  
    # add labels
  if (xlab == paste(col1) & !is.null(axis.col)) {
    xlab <- paste(axis1)
  }
  if (ylab == paste(col2) & !is.null(axis.col)) {
    ylab <- paste(axis2)
  }
  title(xlab = xlab, ylab = ylab, line = 13, cex.axis = 1.2, font.lab = 2)
  
  # x axis names
  xaxis.names <- spec.ax1$Species
  axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
  axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
  if (any(str_count(xaxis.names, " ") == 2))
    axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") == 2), grep("^[A-Z]", xaxis.names))], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[intersect(which(str_count(xaxis.names, " ") == 2), grep("^[A-Z]", xaxis.names))]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[intersect(which(str_count(xaxis.names, " ") == 2), grep("^[A-Z]", xaxis.names))]), sep = "~")), las = 2, cex.axis = 1.05)
  
  # number of specimens in Individual ID
  ax.3 <- table(data$ax1)[match(spec.ax1$Abbreviation, names(table(data$ax1)))]
  ax.3 <- ifelse(is.na(ax.3), 0, ax.3)
  axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, ax.3, cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))
  
  # y axis names
  yaxis.names <- spec.ax2$Species
  axis(2, seq(0,((dim2 - 1)*((1+1/(dim1 - 1)))/dim1), length.out = dim2)[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
  axis(2, seq(0,((dim2 - 1)*((1+1/(dim1 - 1)))/dim1), length.out = dim2)[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
  if (any(str_count(yaxis.names, " ") == 2))
    axis(2, seq(0,((dim2 - 1)*((1+1/(dim1 - 1)))/dim1), length.out = dim2)[intersect(which(str_count(yaxis.names, " ") == 2), grep("^[A-Z]", yaxis.names))], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[intersect(which(str_count(yaxis.names, " ") == 2), grep("^[A-Z]", yaxis.names))]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[intersect(which(str_count(yaxis.names, " ") == 2), grep("^[A-Z]", yaxis.names))]), sep = "~")), las = 2, cex.axis = 1.05)
  
  
  # number of specimens in Definitive ID or count of specimens in axis 2
  if (length(data.long$Person[data.long$Person == data.long$Person[1]]) == nrow(data)) {
    ax.4 <- table(data$ax2)[match(spec.ax2$Abbreviation, names(table(data$ax2)))]
    ax.4 <- ifelse(is.na(ax.4), 0, ax.4)
    axis(4, seq(0,1, length.out = length(yaxis.names)) - 0.002, ax.4, cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))
  } else {
    axis(4, seq(0,((dim2 - 1)*((1+1/(dim1 - 1)))/dim1), length.out = dim2) - 0.002, table(factor(data.long[, col2], levels = spec.ax2$Abbreviation))/length(unique(data.long$Person)), cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))
  }
  
  # make sure the box goes all the way round
  rect(xlim1, ylim2, xlim2, ylim1, col = rgb(0, 0, 0, alpha = 0))
  
  # add grid
  if (grid) {
    abline(v = seq(0,1, length.out = length(xaxis.names))[1:(length(xaxis.names)/grid.sp)]*grid.sp + (1/(length(xaxis.names) - 1))*(grid.sp-0.5), col = grid.col)
    abline(h = seq(0,((dim2 - 1)*((1+1/(dim1 - 1)))/dim1), length.out = dim2)[1:(length(yaxis.names)/grid.sp)]*grid.sp + ((dim2 - 1)*((1+1/(dim1 - 1)))/dim1)/(dim2 - 1)*(grid.sp-0.5), col = grid.col)
  }
  
  if (key) {
    # add key
    par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
    par(mai = c(1, 0.4, 0.8, 0.5))
    names.key <- rep("", 101)
    names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
    key <- rep(1, 101)
    barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", do.call(palette, list(100))), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
    axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
  }
}