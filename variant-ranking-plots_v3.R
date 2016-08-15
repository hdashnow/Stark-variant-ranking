library("ggplot2")
library("reshape2")
library("scales")
library('plyr')
library('dplyr')
library("RColorBrewer")

data_file = "~/Dropbox/Harriet's Documents/MCRI Bioinformatics/variant-ranking-MGHA/data/variantrankingsV3.txt"

csv.colnames = c("Study.ID","HPO.terms","Variant","No.of.variants.available","GPI",
                 "VPI","Combined","Exomiser","CADD","Condel",
                 "CADD.score","Condel.score", 
                 'GERP.rank','GERP.score','PhyloP.rank','PhyloP.score','X','Y','Z')
                  # XYZ are just extra columns containing NAs.

ranking.cols = c("GPI", "VPI", "Combined", "Exomiser","CADD","Condel") # Columns containing ranks
plotting.cols =  c("Combined", "VPI", "Exomiser","CADD","Condel") # Columns to plot
table.cols = c("Combined","Exomiser","CADD","Condel") # Columns for table 1 in paper
my.colours = brewer.pal(6,"Spectral")[c(1,2,3,5,6)] # a colour set that looks okay when printed greyscale

# Swap order of Combined and VPI for the boxplot
boxplot.order = c(2,1,3,4,5)
boxplot.colours = my.colours[boxplot.order]
boxplot.cols = plotting.cols[boxplot.order]

raw.variantrankings = read.delim(data_file, 
                                 comment.char="#", stringsAsFactors=FALSE,
                                 col.names=csv.colnames)

# Remove trailing empty columns
raw.variantrankings[,c('X','Y','Z')] = list(NULL)
# Clean up null values
null_values = c('none', 'None', 'no score')
for (null_val in null_values) {
  raw.variantrankings[raw.variantrankings == null_val] = NA
  }
# Change number columns to numeric type
raw.variantrankings[,4:16] = sapply(raw.variantrankings[,4:16], as.numeric)
# Split out gene names
raw.variantrankings$gene = do.call('rbind', strsplit(as.character(raw.variantrankings$Variant),',',fixed=TRUE))[,1]

data = raw.variantrankings

# Re-rank as number of irrelevent variants looked at before all correct ones found
# -1 to all ranked variants, and -1 for each correct variant ranked higher than it in the same patient

# Re-rank: "adjusted rank is equal to the original rank minus the number of causal variants ranked higher"
re.rank = function(v) {sapply(v, function(x){x - sum(v<x, na.rm = T)})} # Number of non-causal variants above this one in list
data = ddply(data, "Study.ID", function(x) {
  x[,ranking.cols] = apply(x[,ranking.cols], 2, re.rank)
  return(x)
})

# Replace NAs with...
#data[is.na(data)] = 200

#http://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

## Box plot

make.boxplot = function(data, title='') {
  
  pdf('paper/boxplot_v3.pdf', width = 10)
  melted.data = melt(data[,boxplot.cols])
  p = ggplot(melted.data, aes(factor(variable), value, fill=factor(variable))) + 
    geom_boxplot() + scale_y_continuous(trans=reverselog_trans(10), 
                                        breaks = base_breaks()) +
    labs(title=title, x='Ranking method',y="Rank", fill = 'Ranking method') +
    theme_bw() +
    scale_fill_manual(values=boxplot.colours)
  print(p)
  dev.off()
}

title = paste0('Ranks of causal variants, N = ', nrow(data))
make.boxplot(data, title)


## Cumulative density plot

# Calculate my own ecdf so I can control the plotting more finely
my_ecdf = function(x) {
  values = sort(unique(x))
  count = sapply(values, function(y){sum(x == y, na.rm=T)})
  cum.count = cumsum(count)
  return(cbind(values, cum.count))
}

df_ecdf = ddply(melt(data[,plotting.cols]), .(variable), function(x){
  my_ecdf(x$value)
               })

pdf('paper/cumulative_v3.pdf', width = 10)
p = ggplot(df_ecdf,aes(x=values, y=cum.count,colour=variable, linetype=variable)) + 
  geom_step(size = 1.1) +
  scale_x_continuous(trans='log10', breaks = base_breaks()) + 
  labs(x='Rank of the causal variant(s)',y="Number of variants", 
       color = 'Ranking method', 
       linetype = 'Ranking method', 
       title="Number of causal variants found with a given rank")  +
  theme_bw() +
  scale_colour_manual(values=my.colours)
print(p)
dev.off()


## Permutation test

ranks = data[,table.cols]

permute_data = function(df){
  permuted.df = t(apply(df, 1, {function(row) sample(row, length(row), replace = F)}))
  colnames(permuted.df) = names(df)
  return(permuted.df)
}

# Get the ranking method with the lowest mean
lowest_mean = function(df) {
  return(colnames(df)[which.min(apply(df, 2, mean, na.rm  = T))])
}

n_permutations = 999
lowest_mean_results = rep(NA, n_permutations)
for (i in 1:n_permutations) {
  df = permute_data(ranks)
  lowest_mean_results[i] = lowest_mean(df)
}
barplot(table(lowest_mean_results))
lowest_mean(ranks)
p = (sum(lowest_mean_results == "Combined") + 1)/(n_permutations + 1)
p

# Get the ranking method with the lowest mean
topN = function(df, N=5) {
  apply(df, 2, {function(x) sum(x<=N, na.rm  = T)/length(x)})
}

n_permutations = 999
N = 1
topN_results = matrix(NA, nrow = n_permutations, ncol = length(names(ranks)))
colnames(topN_results) = names(ranks)
for (i in 1:n_permutations) {
  df = permute_data(ranks)
  topN_results[i,] = topN(df, N=N)
}
barplot(table(topN_results[,1]))
proportions = topN(ranks, N=N)

calc_p = function(method) {(sum(topN_results[,method] >= topN(ranks)[method]) +1)/(n_permutations + 1)}
p_vals = sapply(colnames(ranks), calc_p)
cbind(proportions, p_vals)[order(proportions, decreasing = T),]


N = 5
topN_results = matrix(NA, nrow = n_permutations, ncol = length(names(ranks)))
colnames(topN_results) = names(ranks)
for (i in 1:n_permutations) {
  df = permute_data(ranks)
  topN_results[i,] = topN(df, N=N)
}
barplot(table(topN_results[,1]))
proportions = topN(ranks, N=N)

calc_p = function(method) {(sum(topN_results[,method] >= topN(ranks)[method]) +1)/(n_permutations + 1)}
p_vals = sapply(colnames(ranks), calc_p)
cbind(proportions, p_vals)[order(proportions, decreasing = T),]


### Get counts, proportions, ranges etc. for various things (for table 1)
# count NAs
apply(data[,table.cols], 2, function(x){ sum(is.na(x)) })
# percent
apply(data[,table.cols], 2, function(x){ sum(is.na(x))/length(x)*100 })

# Mean rank
apply(data[,table.cols], 2, mean, na.rm = T)
# range
apply(data[,table.cols], 2, min, na.rm = T)
apply(data[,table.cols], 2, max, na.rm = T)

# Rank = 1
apply(data[,table.cols], 2, function(x){ sum(x==1, na.rm = T) })
# percent
apply(data[,table.cols], 2, function(x){ sum(x==1, na.rm = T)/length(x)*100 })

# Top 5
apply(data[,table.cols], 2, function(x){ sum(x<=5, na.rm = T) })
# percent
apply(data[,table.cols], 2, function(x){ sum(x<=5, na.rm = T)/length(x)*100 })

# Of the X variants which were not ranked as a result of their inclusion in the gene list-derived GPI, 
# the average rank produced by the VPI was:
VPIs.without.GPIs = data$VPI[is.na(data$GPI)]
length(VPIs.without.GPIs)
mean(VPIs.without.GPIs)

# Condel/CADD scores
apply(data[,c("CADD.score","Condel.score")], 2, mean, na.rm = T)

# Compare combined and VPI along
mean(data$VPI - data$Combined)

# Mean Combined rank for the 39 variants in clinician lists (with GPIs)
Combined.with.GPIs = data$Combined[!is.na(data$GPI)]
mean(Combined.with.GPIs)
# Number of patients with GPIs:
length(unique(data$Study.ID[!is.na(data$GPI)]))

# Total number of patients
length(unique(data$Study.ID))
# Number of patients with multiple variants deemed causative
length(unique(data$Study.ID[duplicated(data$Study.ID)]))

