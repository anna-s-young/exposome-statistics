# Analyzing mixture effects of the untargeted chemical exposome

These snippets of R code are intended to be used within an existing workflow for customized modifications to weighted quantile sum (WQS) regression mixture models. Please find detailed information in our upcoming publication (Young et al, 2025, Submitted to _Exposome_, "Commentary: A Statistical Workflow for Analyzing the Untargeted Chemical Exposome and Metabolome in Epidemiologic Studies Using Mixture Methods").



## Code to Manually Quantize Chemical Exposure Data with Non-Detects

This example code allows you to quantize each exposure variable while separating non-detects (zero values) into their own separate quantile, for a total of 10 deciles. It is modified from the default gwqs_rank() function of the gWQS package in R. It may require modification if your non-detects are indicated in a different way than zeroes or if you want a different type of quantile than deciles.

```
q = 9 # number of non-zero quantiles (zero intensities will be an additional separate one)
data_quantized = data # use data frame where each chemical is in a separate column (observations in rows)
for (chem in chemical_list) { # chemical_list should be a vector of the columns corresponding to exposure variables
  quantiles = unique(quantile(data_quantized[!is.na(data_quantized[,chem]) & data_quantized[,chem]!=0,chem], 
                              probs = seq(0, 1, by=1/q), na.rm=T))
  if (length(quantiles) == 1) {
    quantiles = c(-Inf, quantiles)
  } else {
    quantiles[1] = -Inf
    quantiles[length(quantiles)] = Inf
  }
  binned = cut(data_quantized[!is.na(data_quantized[,chem]) & data_quantized[,chem]!=0,chem], breaks=quantiles, labels=F, include.lowest=T) # q1-->q9
  data_quantized[!is.na(data_quantized[,chem]) & data_quantized[,chem]!=0,chem] = binned
  data_quantized[!is.na(data_quantized[,chem]) & data_quantized[,chem]==0,chem] = 0 # so zeros will be their own quantile of 0
}
# set q=NULL in the gwqs() function so that the WQS model will not automatically quantize your data
```
Here is a (very small) example of the formatting for the data frame called "data" before quantization:

![unquantized_data_example](https://github.com/user-attachments/assets/6098ef00-bd44-4527-93ab-647d17b7718c)

And "data_quantized" after quantization:

![quantized_data_example](https://github.com/user-attachments/assets/309acd10-8f40-428c-ab3d-eca3c1308791)



## Code to Manually Randomize WQS "Repeated Holdouts" by Pairs of Observations

This example code manually defines the WQS repeated holdouts of training versus testing data subsets based on randomized pairs of observations instead of single observations, such that a matched case-control pair of patients (or other pairing type) are always kept within the same data subset.

```
n_RHs = 100 # the number of repeated holdouts you're using
wqs_validation_rows = list()
pairs = unique(data_quantized$PairSet) # use column that indicates the ID of the matched pair
temp = data_quantized
for (r in 1:n_RHs) {
  pairs_partition = sample(c(F,T), size=length(pairs), replace=T, prob=c(1-0.6, 0.6)) # randomly shuffle by pair (~60% validation/testing set)
  which_pairs = pairs[pairs_partition]
  temp$pair_in_partition = F # to reset each iteration
  temp$pair_in_partition[temp$PairSet %in% which_pairs] = T
  wqs_validation_rows[[r]] = temp$pair_in_partition
} 
# input wqs_validation_rows as value for argument "validation_rows" in gwqs() function of gWQS package
```
Example of output (a list of 100 vectors -- one for each repeated holdout -- indicating which samples to include in training vs testing set):
![repeated_holdout_output_example](https://github.com/user-attachments/assets/f202f793-aa18-4090-9758-fb71ca603789)



## Code to Visualize WQS Models with Repeated Holdouts and Random Subsets

This example code can be modified to display results of WQS models using the information on all repeated holdouts of data.

```
## get chemical weights from "model" (gwqs output)

feature_weights = model$final_weights
threshold = 1/nrow(feature_weights)
print(paste("Number of chemicals with average weights above the equi-weight threshold:", nrow(feature_weights %>% filter(Estimate > threshold))))


## load individual repeated holdout data

wmat = data.frame(model$wmat)
wmat$holdout = rownames(wmat)
n_RH = nrow(wmat)
holdouts = wmat %>% pivot_longer(colnames(wmat)[colnames(wmat)!="holdout"], names_to="mix_name", values_to="mean_weight") %>%
  mutate(above_threshold = ifelse(mean_weight > threshold, 1, 0))


## calculate repeated holdout summary stats

holdouts_stats = holdouts %>% group_by(mix_name) %>%
  summarise(num_rh_above_threshold = sum(above_threshold), mean_weight = mean(mean_weight)) %>% ungroup()
holdouts_stats = holdouts_stats %>% mutate(
    percent_rh_above_threshold = num_rh_above_threshold / n_RH * 100,
    contributor = case_when(
      percent_rh_above_threshold >= 90 ~ "probable",
      percent_rh_above_threshold >= 50 ~ "possible",
      percent_rh_above_threshold >= 10 ~ "possibly not",
      percent_rh_above_threshold < 10 ~ "probably not"), 
    mean_weight_above_threshold = ifelse(mean_weight > threshold, "Yes", "No")) %>% 
  arrange(desc(num_rh_above_threshold))

print("Busgang criteria for chemical contributors:")
print(table(holdouts_stats$contributor))


## visualize results

# choose chemicals to display in the graph (depending on how many can fit)
chemicals_to_display = holdouts_stats %>% filter(contributor %in% c("probable")) %>% pull(mix_name)

# set a scale coefficient for the secondary axis in the graph
scale_coeff = max(holdouts_stats %>% filter(mix_name %in% chemicals_to_display) %>% pull(percent_rh_above_threshold)) /
  max(holdouts %>% filter(mix_name %in% chemicals_to_display) %>% pull(mean_weight))

# graph
ggplot(holdouts_stats %>% filter(mix_name %in% chemicals_to_display), aes(y=mix_name, x=percent_rh_above_threshold)) + 
  geom_bar(stat='identity', alpha=0.2, width=0.7, fill="#2c7a83") +
  geom_boxplot(data=holdouts %>% filter(mix_name %in% chemicals_to_display), 
               aes(y=mix_name, x=mean_weight * scale_coeff), outlier.size=0.7, width=0.7, alpha=0.6, color="#2c7a83", fill="#2c7a83") + 
  geom_vline(aes(xintercept=threshold * scale_coeff)) + # equi-weight threshold line
  geom_label(data=data.frame(mean_weight=threshold, mix_name=chemicals_to_display[length(chemicals_to_display)], label="threshold"), aes(x=mean_weight * scale_coeff, y=mix_name, label=label), color="black", fill="white", size=2.3, hjust=0.47, vjust=-0.5) +
  #scale_y_discrete(limits=chemicals_to_display, position="right") + # can add a "labels" parameter with chosen names/annotations of the chemical IDs
  scale_x_continuous(name="% of Holdouts with Weight > Threshold (Barchart)", limits=c(0, max(holdouts_stats %>% filter(mix_name %in% chemicals_to_display) %>% pull(percent_rh_above_threshold))),
                     sec.axis=sec_axis(trans=~./scale_coeff, name="Distribution of Weights (Boxplot)")) +
  theme_minimal() + theme(legend.position="left", panel.grid.minor=element_blank(), axis.text.y=element_text(size=10)) + ylab("Chemical Component") 

```
![example_graph](https://github.com/user-attachments/assets/5d69d171-5e8c-4b6c-87f7-cc9c1a8b11fb)




## References

gWQS R package: https://cran.r-project.org/package=gWQS

Carrico, C., Gennings, C., Wheeler, D. C., & Factor-Litvak, P. (2015). Characterization of Weighted Quantile Sum Regression for Highly Correlated Data in a Risk Analysis Setting. Journal of Agricultural, Biological, and Environmental Statistics, 20(1), 100–120. https://doi.org/10.1007/s13253-014-0180-3

Tanner, E. M., Bornehag, C.-G., & Gennings, C. (2019). Repeated holdout validation for weighted quantile sum regression. MethodsX, 6, 2855–2860. https://doi.org/10.1016/j.mex.2019.11.008

Curtin, P., Kellogg, J., Cech, N., & Gennings, C. (2021). A random subset implementation of weighted quantile sum (WQSRS) regression for analysis of high-dimensional mixtures. Communications in Statistics - Simulation and Computation, 50(4), 1119–1134. https://doi.org/10.1080/03610918.2019.1577971


