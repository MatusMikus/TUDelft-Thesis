library(ggplot2)
library(ggallin)
library(scales)
library(grid)
library(data.table)
library(dplyr)
library(gdata)
#library(ggbreak)

setwd("C:/Users/matus/Desktop/TUDelft/Cluster Backup/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenVisR")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  print(list(...))
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

make_plot = function(pipeline, color, title,  filter_outliers=TRUE) {
  outliers = boxplot(pipeline$LENGTH, plot = FALSE)$out #naive implementation but enough for our purposes
  print(outliers)
  new_ppln = pipeline #in case we don't want to filter outliers
  if(filter_outliers){
    print("filtering outliers")
    new_ppln = pipeline[!(pipeline$LENGTH %in% outliers), ]    
  }
  df = data.frame(length = new_ppln$LENGTH, count = new_ppln$COUNT)
  plt = ggplot(df,aes(x=length, y=count)) + 
    ggtitle(title) +
    geom_bar(aes(x=length, y=count), stat="identity",colour=color, lwd=0.5) +
    scale_y_continuous(breaks=c(0,5,20,100,250,500,1000,5000,20000), trans=pseudolog10_trans, expand=c(0,0), name="variant count") + 
    scale_x_continuous(name="variant length (bp)")
  print(plt)
  return(plt)
}

adjust_plots = function(..., plotlist=NULL, file, cols=1, layout=NULL, maxValue = 50000){
  plots = c(list(...))
  newPlots = c()
  min = 0
  max = 0
  for (plot in plots) {
    lims = layer_scales(plot)$x$get_limits()
    min = min(min,lims[1])
    max = max(max,lims[2])
  }
  min = -maxValue
  max = maxValue

  for (i in 1:length(plots)) {
    plots[[i]] = plots[[i]] + 
      xlim(min,max) + 
      scale_x_continuous(trans=pseudolog10_trans,breaks=c(-10000,-2000,-400,-100,100,400,2000,10000),limits=c(min,max),name="variant length (bp)")
  }

  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
  
}


#adjust_plots(plt_mhc_dv, plt_mhc_pbsv, plt_mhc_sniffles, maxValue = 10000)

deepvariant_files = list.files(path = "lengths/deepvariant2/",
                               recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
deepvariant_array <- rbindlist(sapply(deepvariant_files, fread, simplify = FALSE),
                use.names = TRUE, idcol = "FileName")

pbsv_files = list.files(path = "lengths/pbsv/",
                               recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
pbsv_array <- rbindlist(sapply(pbsv_files, fread, simplify = FALSE),
                               use.names = TRUE, idcol = "FileName")

#sniffles_files = list.files(path = "lengths/sniffles/",
#                        recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
#sniffles_array <- rbindlist(sapply(sniffles_files, fread, simplify = FALSE),
#                        use.names = TRUE, idcol = "FileName")


sniffles_len_files = list.files(path = "sniffles_lengths/sniffles",
                                recursive = TRUE, pattern = "*.txt", full.names = TRUE)

sniffles_len_array <- rbindlist(sapply(sniffles_len_files, fread, simplify = FALSE),
                                use.names = TRUE, idcol = "FileName")


lengths_sniffles_all = sniffles_len_array %>% group_by(V1) %>% summarise(COUNT = n())
colnames(lengths_sniffles_all) = c("LENGTH","COUNT")

plt_sniffles = make_plot(lengths_sniffles_all, "red", "Sniffles")

#add together all the ones with the same length

lengths_dv = deepvariant_array %>% group_by(LENGTH) %>% filter(abs(LENGTH) >= 10) %>% summarise(COUNT = sum(COUNT)) 
lengths_pbsv = pbsv_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))
lengths_sniffles = sniffles_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))

#barplot(avgmat, log="y", col=c("green","red", "blue"), beside=TRUE, legend = c("DeepVariant", "pbsv", "sniffles"))

plt_dv = make_plot(lengths_dv, "blue", "A")
plt_pbsv = make_plot(lengths_pbsv, "darkgreen", "B")
plt_sniffles = make_plot(lengths_sniffles_all, "red", "C")

adjust_plots(plt_dv, plt_pbsv, plt_sniffles, maxValue = 10000)
adjust_plots(plt_dv, plt_pbsv, plt_sniffles, maxValue = 1000)

#Make a new plot with only short lengths (500)

lengths_dv_short = deepvariant_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))
lengths_pbsv_short = pbsv_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))
lengths_sniffles_short = sniffles_len_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))

plt_dv_short = make_plot(lengths_dv_short, "blue", "DeepVariant")
plt_pbsv_short = make_plot(lengths_pbsv_short, "darkgreen", "PBSV")
plt_sniffles_short = make_plot(lengths_sniffles_short, "red", "Sniffles")

adjust_plots(plt_dv_short, plt_pbsv_short, plt_sniffles_short, maxValue = 410)

#
#bar chart for lengths
#

vector_seq = seq(from = -1000, to = 1000 , by = 100)
bins_dv = lengths_dv %>% mutate(length_bin = cut(LENGTH, breaks=vector_seq))
bins_g_dv = tapply(bins_dv$COUNT, bins_dv$length_bin, sum)
barplot(height = bins_g_dv)

bins_pbsv = lengths_pbsv %>% mutate(length_bin = cut(LENGTH, breaks=vector_seq))
bins_g_pbsv = tapply(bins_pbsv$COUNT, bins_pbsv$length_bin, sum)
barplot(height = bins_g_pbsv)

bins_sniffles = lengths_sniffles_all %>% mutate(length_bin = cut(LENGTH, breaks=vector_seq))
bins_g_sniffles = tapply(bins_sniffles$COUNT, bins_sniffles$length_bin, sum)
barplot(height = bins_g_sniffles)

combined_2 = t(data.frame(bins_g_dv, bins_g_pbsv,bins_g_sniffles))

barplot(as.matrix(combined_2), beside=TRUE)



#create a bar chart for regions:
dat <- read.table('region_stats/deepvariant2/HG002.txt', skip = 1, stringsAsFactors = FALSE)
dat_h = dat[seq(1, nrow(dat), 2),]
dat_d = dat[seq(2, nrow(dat), 2),]
#dat_d$names = dat_h
df = data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) = dat_h
df[nrow(df) + 1,] = dat_d

extract_region_values = function(file){
  dat <- read.table(file, skip = 1, stringsAsFactors = FALSE)
  dat_h = dat[seq(1, nrow(dat), 2),]
  dat_d = dat[seq(2, nrow(dat), 2),]
  dat_d$names = dat_h
  df = data.frame(matrix(ncol = 5, nrow = 0)) #create empty data frame
  colnames(df) = dat_h # add column names from extracted headers
  df[nrow(df) + 1,] = dat_d #add new row
  return(df)
}

deepvariant_regions = list.files(path = "region_stats/deepvariant2/",
                               recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
deepvariant_r_array <- rbindlist(sapply(deepvariant_regions, extract_region_values, simplify = FALSE),
                               use.names = FALSE)

pbsv_regions = list.files(path = "region_stats/pbsv/",
                                 recursive = TRUE, pattern = "*.txt", full.names = TRUE)

pbsv_r_array <- rbindlist(sapply(pbsv_regions, extract_region_values, simplify = FALSE),
                                 use.names = FALSE)

sniffles_regions = list.files(path = "region_stats/sniffles/",
                          recursive = TRUE, pattern = "*.txt", full.names = TRUE)

sniffles_r_array <- rbindlist(sapply(sniffles_regions, extract_region_values, simplify = FALSE),
                          use.names = FALSE)

#Second plot 

combined = gdata::combine(deepvariant_r_array, pbsv_r_array,sniffles_r_array)
combined = as.data.frame(do.call(cbind,(combined)))
combined[] <- lapply(combined, function(x) as.numeric(as.character(x)))

names = colnames(combined) 
rows = combined[!colnames(combined) %in% c("source")]

avg_combined = combined %>%
  group_by(source) %>%
  summarise_all(list(mean))

avg_combined = avg_combined[!colnames(avg_combined) %in% c("source")]

avg_mat = matrix(avg_combined)
avg_mat = mapply(avg_mat, FUN=as.numeric)
avgmat = avg_mat+1
colnames(avgmat) = colnames(avg_combined)


#move legend
barplot(avgmat, log="y", ylab = "Average variants per genome", xlab="Region",
        col=c("green","red", "blue"), beside=TRUE,
        legend = c("DeepVariant", "PBSV", "Sniffles"))


#create a bar chart for DP:

extract_dp_values = function(file){
  dat_dp =  read.table(file, skip = 0, stringsAsFactors = FALSE)
  dat_dp_h = dat_dp[seq(1, nrow(dat_dp), 2),]
  dat_dp_d = dat_dp[seq(0, nrow(dat_dp), 2),]
  dat_dp_d$names = dat_dp_h
  df = data.frame(matrix(ncol = 7, nrow = 0))
  colnames(df) = dat_dp_h
  df[nrow(df) + 1,] = dat_dp_d
  return(df)
}


deepvariant_dps = list.files(path = "dp_stats/deepvariant2/",
                                 recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
deepvariant_dp_array <- rbindlist(sapply(deepvariant_dps, extract_dp_values, simplify = FALSE),
                                 use.names = FALSE)

pbsv_dps = list.files(path = "dp_stats/pbsv/",
                          recursive = TRUE, pattern = "*.txt", full.names = TRUE)

pbsv_dp_array <- rbindlist(sapply(pbsv_dps, extract_dp_values, simplify = FALSE),
                          use.names = FALSE)

#Second plot 

colnames = c(colnames(deepvariant_dp_array), "source")
combined_dp = gdata::combine(deepvariant_dp_array, pbsv_dp_array)


colnames(combined_dp) = colnames
combined_dp = as.data.frame(do.call(cbind,(combined_dp)))
combined_dp[] <- lapply(combined_dp, function(x) as.numeric(as.character(x)))

names_dp = colnames(combined_dp) 
rows_dp = combined_dp[!colnames(combined_dp) %in% c("source")]

avg_combined_dp = combined_dp %>%
  group_by(source) %>%
  summarise_all(list(mean))

avg_combined_dp = avg_combined_dp[!colnames(avg_combined_dp) %in% c("source")]
names=c("10-","10-19","20-29","30-49","50-99","100-199","200+")
colnames(avg_combined_dp) = names
avg_dp_mat = matrix(avg_combined_dp)
avg_dp_mat = mapply(avg_dp_mat, FUN=as.numeric)
avgmat = avg_mat+1
colnames(avg_dp_mat) = names

barplot(avg_dp_mat, log="y", 
        ylab="Average variants per genome", xlab="Read Depth of variant", 
        col=c("green","red"), beside=TRUE)
legend("topright",
       pch = 15,
       legend = c("DeepVariant", "pbsv"), 
       col = c("green","red")
     )


#Investigation:

sniffles_regions = list.files(path = "region_stats/sniffles/",
                              recursive = TRUE, pattern = "*.txt", full.names = TRUE)

sniffles_r_array <- rbindlist(sapply(sniffles_regions, extract_region_values, simplify = FALSE),
                              use.names = FALSE)




##### MHC

pbsv_mhc_array = read.delim("mhc/pbsv_region_lengths.txt")
lengths_mhc_pbsv = pbsv_mhc_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))

dv_mhc_array = read.delim("mhc/deepvariant_region_lengths.txt")
lengths_mhc_dv = dv_mhc_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))

sniffles_mhc_array = read.delim("mhc/sniffles_fixed_lengths.txt")
sniffles_mhc_lengths = sniffles_mhc_array %>% group_by(X38) %>% summarise(COUNT = n())
colnames(sniffles_mhc_lengths) = c("LENGTH","COUNT")
sniffles_mhc_lengths = sniffles_mhc_lengths %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))

plt_mhc_dv = make_plot(lengths_mhc_dv, "blue", "A", filter_outliers = FALSE)
plt_mhc_pbsv = make_plot(lengths_mhc_pbsv, "darkgreen", "B",filter_outliers = FALSE)
plt_mhc_sniffles = make_plot(sniffles_mhc_lengths, "red", "C", filter_outliers = FALSE)

adjust_plots(plt_mhc_dv, plt_mhc_pbsv, plt_mhc_sniffles, maxValue = 10000)



deepvariant_reg_files = list.files(path = "lengths/deepvariant2/",
                               recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
deepvariant_reg_array <- rbindlist(sapply(deepvariant_files, fread, simplify = FALSE),
                               use.names = TRUE, idcol = "FileName")

pbsv_reg_array = list.files(path = "lengths/pbsv/",
                        recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
pbsv_reg_array <- rbindlist(sapply(pbsv_files, fread, simplify = FALSE),
                        use.names = TRUE, idcol = "FileName")

sniffles_files = list.files(path = "lengths/sniffles/",
                        recursive = TRUE, pattern = "*.txt", full.names = TRUE)

# Read all the files and create a FileName column to store filenames
#sniffles_array <- rbindlist(sapply(sniffles_files, fread, simplify = FALSE),
#                        use.names = TRUE, idcol = "FileName")


sniffles_len_files = list.files(path = "sniffles_lengths/sniffles",
                                recursive = TRUE, pattern = "*.txt", full.names = TRUE)

sniffles_len_array <- rbindlist(sapply(sniffles_len_files, fread, simplify = FALSE),
                                use.names = TRUE, idcol = "FileName")


lengths_sniffles_all = sniffles_len_array %>% group_by(V1) %>% summarise(COUNT = n())
colnames(lengths_sniffles_all) = c("LENGTH","COUNT")

plt_sniffles = make_plot(lengths_sniffles_all, "red", "Sniffles")

#add together all the ones with the same length

lengths_dv = deepvariant_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))
lengths_pbsv = pbsv_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))
lengths_sniffles = sniffles_array %>% group_by(LENGTH) %>% summarise(COUNT = sum(COUNT))

plt_dv = make_plot(lengths_dv, "blue", "DeepVariant")
plt_pbsv = make_plot(lengths_pbsv, "darkgreen", "PBSV")
plt_sniffles = make_plot(lengths_sniffles_all, "red", "Sniffles")

adjust_plots(plt_dv, plt_pbsv, plt_sniffles, maxValue = 10000)
adjust_plots(plt_dv, plt_pbsv, plt_sniffles, maxValue = 1000)

#Make a new plot with only short lengths (500)

lengths_dv_short = deepvariant_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))
lengths_pbsv_short = pbsv_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))
lengths_sniffles_short = sniffles_len_array %>% group_by(LENGTH) %>% filter(LENGTH <= 500 && LENGTH >= -500) %>% summarise(COUNT = sum(COUNT))

plt_dv_short = make_plot(lengths_dv_short, "blue", "DeepVariant")
plt_pbsv_short = make_plot(lengths_pbsv_short, "darkgreen", "PBSV")
plt_sniffles_short = make_plot(lengths_sniffles_short, "red", "Sniffles")

adjust_plots(plt_dv_short, plt_pbsv_short, plt_sniffles_short, maxValue = 410)









