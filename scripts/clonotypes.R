#!/usr/bin/Rscript

# Script for the haplotype clustering image
usePackage <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) {
    return (TRUE)
  }

  install.packages(package, repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

usePackage.bio <- function(package){

  library.path <- .libPaths()[1];

  if(!is.element(package, installed.packages(lib.loc=library.path)[,1]))
    BiocManager::install(package, lib=library.path, update=FALSE)

  # BiocManager::install(package, update=FALSE)
  # install.packages(package, repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

# For command line parsing
usePackage("argparse");

# create parser object
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default=NULL, help="input distribution matrix file"); # Distribution matrix file
parser$add_argument("-o", "--output", default=NULL, type="character");

opt <- parser$parse_args();

# Check options
if (is.null(opt$input) || is.null(opt$output)){
  print_help(parser)
  stop("Missing argument(s)!", call.=FALSE)
}

# Install packages
usePackage("ggplot2");
usePackage("ggdendro"); # To plot dendrograms 
usePackage("tidyr");
usePackage('gridExtra');
usePackage("dplyr");
usePackage("cowplot");  # For plot_grid
usePackage("RColorBrewer");
usePackage("BiocManager");
usePackage("showtext");

options(dplyr.print_max = 1e9)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  # cormat <-cormat[hc$order, hc$order]
  nameList <- rownames(cormat);
  nameList <- nameList[hc$order];
  # cormat <- nameList;
  return(nameList);
}

# Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

ggplot.font <- '';
border.color <- 'black';
border.size <- 0.30;
font_add(family = "Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()

# Load dataset and rename columns
data.frequency <- read.csv(opt$input, header = TRUE, sep="\t");
names(data.frequency) <- c(
  'id',
  'locus',
  'cdr3',
  'cdr3_aa',
  'v_family',
  'd_family',
  'j_family',
  'cdr3_length',
  'type'
);

# Calculate frequency
data.frequency <- data.frequency %>%
  mutate(label = paste(v_family, j_family, sep = ";")) %>%
  group_by(label, type) %>%
  summarise(count = n());

# Append missing elements to fill the tables
typeList <- unique(data.frequency %>% distinct(type) %>% pull());
labelList <- sort(unique(data.frequency %>% distinct(label) %>% pull()));
data.frequency.empty <- expand.grid(type=typeList, count=NA, label=labelList);
data.frequency.missing <- anti_join(data.frequency.empty, data.frequency, by=c("type", 'label'));
if(nrow(data.frequency.missing) > 0){
  data.frequency.missing$count <- 0;
  data.frequency <- full_join(data.frequency, data.frequency.missing, by=c("type", 'label', 'count'));
}

# Calculate frequency
data.frequency <- data.frequency %>%
  group_by(type) %>%
  mutate(total=sum(count), freq = count / total * 100) %>%
  mutate(freq = ifelse(count==0, NA, freq)) %>%
  ungroup();

# Sort data by label alphabetical order
data.frequency$label <- factor(data.frequency$label, levels=sort(unique(data.frequency$label)));

# Generate a correlation plot for each group
groupList <- data.frequency %>% distinct(type) %>% pull();
group.correlation <- matrix(data=NA,nrow=length(groupList), ncol=length(groupList));
colnames(group.correlation) <- groupList;
rownames(group.correlation) <- groupList;
correlationTest <- c();
for(i in 1:length(groupList)){

  # Generate the vector for A
  groupA <- data.frequency %>% filter(type == groupList[i]) %>% arrange(label) %>% select(freq) %>% pull() %>% replace_na(0);
  for(j in i:length(groupList)){

    # Generate the vector for B
    groupB <- data.frequency %>% filter(type == groupList[j]) %>% arrange(label) %>% select(freq) %>% pull() %>% replace_na(0);

    # Calculate pearson correlation
    correlation <- cor(groupA, groupB, method = "pearson");
    correlation.test <- cor.test(groupA, groupB, method = "pearson");
    r2 <- cor.test(groupA, groupB, method = "pearson")$p.value
    r3 <- symnum(r2, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*",""));
    correlationTest <- c(correlationTest, r3);

    cat(paste('Correlation between', groupList[i], 'and', groupList[j], ':', correlation, '\n'));
    if(is.na(correlation)){
      correlation <- 0;
    }
    group.correlation[i,j] = correlation;
    group.correlation[j,i] = correlation;
  }
}

# Plot the haplotype
p.heatmap <- ggplot(data.frequency) + 
  geom_tile(mapping=aes(x=label, y=type, fill =type, alpha=freq), colour=border.color, size=border.size) +
  facet_grid(rows = vars(type), switch="both", scales = "free", space = "free") +
  scale_fill_manual(labels = c("Sn 2nd dose", "Sn 3rd dose", 'Seropositives'), values = c("#C8DEF9", "#004080", '#A00000')) +
  guides(colour='none', fill=guide_legend(nrow = 1, title.position = "left")) +
  xlab('Clonotypes') + 
  theme_void() +
  theme(
        text = element_text(family = ggplot.font),
        axis.text.x=element_text(angle = 90, size=4, hjust = 1, margin = margin(t=0, r=0, b=0, l=0)),
        axis.text.y=element_blank(),
        strip.text.y.left=element_text(angle = 0),
        strip.text.y=element_text(angle = 0, size=4, hjust = 0.5, margin = margin(t=0, r=0, b=0, l=0)),
        panel.spacing.y=unit(0, "line"),
        legend.key.size = unit(4, 'mm'),
        strip.placement = "outside",
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(l=0, t=0, b=0, r=0, 'line'),
        axis.line.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank()
        );

# Get cluster legend
p.heatmap.legend  <- get_legend(p.heatmap + theme(
  legend.title=element_blank(), 
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box.just = "bottom", 
  legend.text=element_text(size=8)));

# Remove the legend before plotting
p.heatmap <- p.heatmap + theme(legend.position='none');

# Open a PDF for plotting; units are inches by default
cat("Rendering data..", "\n");

# Plot distinct motifs bar plot 
# This graph shows the number of distinct peptide motives recognised by the same haplotype
# Shared motives appears multiple times in the haplotype, while unique only one time from one specific allele
# allele_scaling <- 8;
p.stacked <- ggplot(data.frequency, aes(x=label, y=freq, fill=type)) +
  geom_bar(stat='identity', colour=border.color, size=border.size) +  
  ylab('Frequency (%)') + 
  scale_fill_manual(labels = c("Sn 2nd dose", "Sn 3rd dose", 'Seropositives'), values = c("#C8DEF9", "#004080", '#A00000')) +
  theme_bw() + 
  theme(
    text = element_text(family = ggplot.font),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_blank(),
    axis.text.y = element_text(size=8),
    axis.line=element_line(color='black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    legend.direction="horizontal",
    legend.title=element_blank(),
    legend.key.size = unit(4, 'mm'),

    plot.margin = unit(c(0, 0, 0, 0), "cm")
  );

# Get legends
p.stacked.legend <- get_legend(p.stacked +
  theme(
    legend.title=element_blank(), 
    legend.direction = "horizontal",
    legend.justification="center",
    legend.box.just = "bottom", 
    legend.text=element_text(size=8),
    plot.margin = margin(t=0, r=0, b=0, l=0, unit="mm")
  ));

# Remove the legend before plotting
p.stacked <- p.stacked + theme(legend.position='none');

# Reorder the correlation matrix
group.correlation.order <- reorder_cormat(group.correlation);

# group.correlation <- get_upper_tri(group.correlation);
# print(group.correlation);
group.correlation <- group.correlation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  gather("Var2", "value", -Var1);# %>%

group.correlation$Var1 <- factor(group.correlation$Var1, levels = group.correlation.order);
group.correlation$Var2 <- factor(group.correlation$Var2, levels = group.correlation.order);

# Create an heatmap of group correlation
p.correlation <- ggplot(group.correlation, aes(Var1, Var2, fill = value))+
  geom_tile(color = border.color, linetype = 1, size=border.size) +
  geom_text(aes(label = round(value, 2)), size=3) +
  scale_fill_gradient2(
    low = "blue",
    high = "red", 
    mid = "white", 
    midpoint = 0.5, 
    limit = c(0,1), 
    space = "Lab",
    na.value = "white",
  ) +
  ylab('Correlation') + 
  theme_void() +
  theme(
    text = element_text(family = ggplot.font),
    legend.title=element_blank(), 
    axis.text.x=element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
    axis.text.y=element_text(size = 8),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.border=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=2, r=2, b=0, l=2, unit="mm"),
    axis.title.y = element_text(angle = 90, vjust = 1, size = 8, hjust = 0.5),
    legend.text=element_text(size=6),
    legend.key.size = unit(14, 'pt'),
    legend.justification = c(0, 1)
  ) +
  coord_fixed();


grid <- plot_grid(
    plot_grid(
        p.stacked,
        p.heatmap,
        nrow=2,
        ncol=1,
        align = 'v', axis = 'l',
        rel_heights=c(1.0, 0.8)
    ),
    plot_grid(
      p.heatmap.legend,
      p.correlation,
      ncol=2
    ),
  
  nrow=2,
  rel_heights=c(1.0,0.5)
) + theme(plot.margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="cm"));

# Save the final plot
ggsave(filename=opt$output, plot=grid, units="mm", width=170, height=120, dpi=300);
showtext_auto(FALSE)

# Often the most useful way to look at many warnings:
summary(warnings())
