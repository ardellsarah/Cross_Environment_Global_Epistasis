# DESCRIPTION:  This Script takes the finalized mutation effect data and creates all main text and supplementary figures
# Other than S02 B -> this is created in the script that processes raw barcode counts to fitness effects


# PACKAGES USE 
library(ggplot2)
library(ggpubr)
library('MASS')
library(viridis)
library(dplyr)
library( tidyr )
library(stats4) 
library(tidyverse)
library(lmtest)
library(moments)
library(purrr)
library(broom)
library(ggridges)
library(ggeffects)
library(microseq) # use to get sequences out from fasta files on computer
library(stringdist) # calc string distances 
library(R.utils) # unzip fastq.gz files
library(bioseq)
library(readxl)


## ALL FILES NEEDED in same folder 


##  Where to save files and core files 
setwd('~/Desktop/TnSeq_Code/') # alter to be the name of the folder with all the data in it
figSave = 'Figures_Revision/' # !! Create this folder in the above directory
fileSave = 'Files_Revision/'  # !! Create this folder in the above directory

neutralMutIDs = c(102,51,99,91,6)
minNumStrains = 4 # min num dat points to do a regression on 


## All Data Need 
df_s_ests_wGR = read.csv( 'df_s_ests_wGR_revision.csv')
df_s_ests_wGR = df_s_ests_wGR[,-1]
growthRate_data = read.csv('growthRate_data_revision.csv')
growthRate_data = growthRate_data[,-1]
JohnsonDat = read.csv(paste0( 'JohnsonDFEmeanDat.csv'))
Mut_full_info_Use = read.csv( 'Mut_full_info_Use.csv')
Mut_full_info_Use = Mut_full_info_Use[, !names(Mut_full_info_Use) %in% 'X']
# df_allSests_allBCs = read.csv('allBC_s_ests.csv')
bloom_temp_pH_qtls = read.csv('Bloom_QTLS_tempPH.csv')
allGeneDiffs_wide = read.csv('elife-27167-supp1-v2.csv')
allExpGR_data = read.csv('allExpGR_data.csv')
Milo_S_est_wide = read.csv('Milo_S_est.csv')
milo_slopes_ints = read.csv('milo_slopes_ints.csv')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####



######## Figure 1 : Beneficial and Deleterious Mutations ##### 

# Prep data: only keep s-ests for mutations with sufficient data
df_mutsCanFitLine = df_s_ests_wGR %>%
  group_by(env, Mut_ID) %>%
  mutate(totStrainsIn = length(Strain))

df_mutsCanFitLine = df_mutsCanFitLine %>%
  group_by(Mut_ID)%>%
  mutate(inEnoughEnvs = sum(unique(totStrainsIn)>=minNumStrains))

df_mutsCanFitLine = subset(df_mutsCanFitLine, totStrainsIn>=4 & !(Mut_ID %in% neutralMutIDs) &inEnoughEnvs>1)
# remove neutral muts because want to focus on non-putative neutral, and remove those that dont have enough strains in at least 2 envs (so can actually comp env dependent slope)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

## ~~ Plot: Prop Ben and Delt Muts vs GR ####

prop_ben_delt = df_mutsCanFitLine %>%
  group_by(Strain, env, Mean_GR) %>%
  summarize(PropBen = sum(mutCall == 'Beneficial')/length(mutCall),
            PropDelt = sum(mutCall == 'Deleterious')/length(mutCall))

prop_ben_delt_g = gather(prop_ben_delt, key = 'type', value = 'Proportion',PropBen:PropDelt )


plotPben_Delt = function(envID) 
{
  colorDF = data.frame(env = c("30SC3" ,"30SC5", "30SC7" ,"37SC3", "37SC5" ,"37SC7" ,"YPD"  ), color = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'))
  
  colorChoice = subset(colorDF, env == envID)$color
  
  data = subset(prop_ben_delt_g, env == envID)
  data$type = factor(data$type, levels = c('PropBen', 'PropDelt'))
  return(  ggplot(data,aes(x = Mean_GR,y = Proportion, color = type, fill = type, linetype = type))+
             geom_point(aes(shape = type), alpha = 0.6, size  = 1.2)+
             geom_smooth(method = 'lm', formula = 'y~x', lineend = 'round', se = F)+
             scale_color_manual(values = c(colorChoice,colorChoice))+
             scale_fill_manual(values = c(colorChoice, 'white'))+
             scale_shape_manual(values = c(17,25))+
             scale_linetype_manual(values = c('solid', 'dashed'))+
             theme_classic() +
             xlab('Growth Rate (1/hr)')+
             ylab('Proportion in DFE')+
             #ggtitle(paste0(envID))+
             scale_y_continuous(expand = c(0,0), limits = c(-0.2,0.68), breaks = c(0, 0.3, 0.6))+
             scale_x_continuous(  breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4))+
             coord_cartesian( ylim=c(0,0.68)) +
             theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
                   #axis.ticks.length=unit(.04, "cm"),
                   axis.title = element_blank(),
                   axis.text.x = element_text(size=8,color = 'black', family = 'Helvetica'),
                   axis.text.y = element_blank(),
                   plot.title = element_text(size=8,color = 'black', family = 'Helvetica'),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) )
  
}



allPben_delt_plots = map(unique(prop_ben_delt_g$env), ~plotPben_Delt( .x))

allPben_delt_plots_a = ggarrange(plotlist = allPben_delt_plots, nrow = 2, ncol = 3)

#ggsave(paste0(figSave,'allPben_delt_plots_a.pdf' ),allPben_delt_plots_a, width = 11.2, height = 7.55, unit = 'cm')






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure 2 : Stats of Slopes and Intercepts (and relevant supplement) ####


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
## ~~ Analysis: Fit variable slopes model ####

# fit eq (1) for each mutation (e.g, get 1 slope and 1 int per env, save Rsq/pvalues)
perMut_6Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR * env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model (eq (1))
perMut_6Slope6int_fits = perMut_6Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_6Slope6int_fits)[2] = paste0(names(perMut_6Slope6int_fits)[2], '_6s6i')
perMut_6Slope6int_Rsq_p = perMut_6Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_6Slope6int_Rsq_p)[2:5] = paste0(names(perMut_6Slope6int_Rsq_p)[2:5], '_6s6i')


# also fits eq (1) for each mutation, but saves Rsq/p for the regression in each environment (or the slopes/ints)
perMut_Env_6Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -c(Mut_ID ,env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR , .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance) 

perMut_Env_6Slope6int_slopes_ints = df_mutsCanFitLine %>%
  nest(data = -c(Mut_ID ,env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR , .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) 
# this gives all slopes and ints 
sub = perMut_Env_6Slope6int_slopes_ints[, c('env', 'Mut_ID', 'term', 'estimate')]
perMut_Env_6Slope6int_slopes_ints_spread = spread(sub, term, estimate)
names(perMut_Env_6Slope6int_slopes_ints_spread) = c('env', 'Mut_ID', 'intercept','slope')
remove(sub)


## ~~ Plot: Histogram of Slopes and Intercepts ####

# TO make y-axis fraction, just specify bindiwth and multiple y value by binwidth 
hist_mutSlopes = ggplot(data = perMut_Env_6Slope6int_slopes_ints_spread, aes(x = slope,color = env)) + 
  stat_binline(binwidth = 0.12, aes(color = env,y = after_stat(density)*0.12),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.5, lineend = 'round')+
  theme_classic()+
  xlab('Slope')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')
#ggsave(paste0(figSave,'hist_mutSlopes.pdf' ),hist_mutSlopes, width = 5, height = 5.4, unit = 'cm')




hist_mutInts = ggplot(data = perMut_Env_6Slope6int_slopes_ints_spread, aes(x = intercept,color = env)) + 
  stat_binline(binwidth = 0.03, aes(color = env,y = after_stat(density)*0.03),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.55, lineend = 'round')+
  theme_classic()+
  xlab('Intercept')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')

#ggsave(paste0(figSave,'hist_mutInts.pdf' ),hist_mutInts, width = 5, height = 5.4, unit = 'cm')



# save parameters of slope and int distributions in each env
# allSlopeIntDistParams = perMut_Env_6Slope6int_slopes_ints_spread %>% 
#   group_by(env) %>%
#   summarize(SlopeDistMean = mean(slope, na.rm = T), SlopeDistSD = sd(slope, na.rm = T), , SlopeDistSkew = skewness(slope),
#             intDistMean = mean(intercept, na.rm = T), intDistSD = sd(intercept, na.rm = T), , intDistSkew = skewness(intercept))

#write.csv(allSlopeIntDistParams, 'allSlopeIntDistParams.csv')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

## ~~ Plot:  Slope Int Corr ####
slopeIntCorr = ggplot(perMut_Env_6Slope6int_slopes_ints_spread, aes(x = slope, y = intercept, color = env, group = env))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_point( alpha = 0.5, size = 2, shape = 16)+
  geom_smooth(method = 'lm',formula = 'y~0+x', se = F, lineend = 'round')+
  xlab('Slope') +  ylab('Intercept')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave,'slopeIntCorr.pdf' ),slopeIntCorr, width = 5, height = 5.4, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure 3 : Microscopic Epistasis Regressions - all Mutations ####

## first get mutation name info 
df_mutsCanFitLine$Mut_ID = as.numeric(df_mutsCanFitLine$Mut_ID)
df_mutsCanFitLine = inner_join(df_mutsCanFitLine,Mut_full_info_Use, by = 'Mut_ID' )

pSigDiffBins = data.frame(min = c(-1,0,0.2,0.4), max = c(0,0.2,0.4,0.6), bin = c(1,2,3,4))
yaxisScales = data.frame(mMin = c(1,57,79),mMax =  c(56,78,94), scaleMin = c(-0.09,-0.09,-0.09), scaleMax = c(0.09,0.09,0.09))
  # here, just using same scales for y-axis for all plots 
plotExample_includepairwiseSlopeDiff = function(rowID) 
{
  m = as.numeric(joinedRanks[rowID,'Mut_ID']) # get the mutation ID assocaited with the row of interest
  
  ## custom y axis scale for plot
  yaxisScales$isIn = rowID >= yaxisScales$mMin & rowID <= yaxisScales$mMax
  ymin = yaxisScales$scaleMin[yaxisScales$isIn]
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  fitted = subset(perMut_6Slope6int_fits, Mut_ID == m)$fit_6s6i[[1]] # get the variable slope model regression
  data = subset(df_mutsCanFitLine, Mut_ID == m)
  data$predS = fitted$fitted.values # note this only works because data here is subsetted from exact same DF in exact same way as when do the LM
  data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  mutName = subset(Mut_full_info_Use, Mut_ID ==m)$Gene.Use
  
  colorScale_pSigDiff = rev(c('black', 'black', 'black', 'black'))
  
  propSigDiff = subset(allPairwiseSlopeDiffs, Mut_ID == m)$PropSigDiff
  pSigDiffBins$isIn = propSigDiff> pSigDiffBins$min & propSigDiff<=pSigDiffBins$max
  
  colorChoice = colorScale_pSigDiff[pSigDiffBins$isIn]
  
  
  return(  ggplot(data,aes(x = Mean_GR))+
             geom_hline(yintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5/2)+
             # geom_point(aes(y = avgS, color = env), alpha = 0.2, size  = 0.4)+
             geom_line(aes(y = predS, color = env),   linewidth = 0.5)+
             scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = F)+
             theme_void() +
             xlab('Growth Rate (1/hr)')+
             ylab('Selection Coefficient (1/hr)')+
             # ggtitle(paste0(mutName))+
             annotate(geom="text", x=-Inf, y=-Inf, label=paste0(mutName),
                      color=colorChoice, size = 1.945, vjust = -0.6, hjust = -0.1)+
             #scale_y_continuous( breaks = c(-0.12,-0.09, -0.06,-0.03, -0.01, 0.00, 0.01,0.03, 0.06, 0.09, 0.12))+
             scale_y_continuous(limits = c(ymin,ymax), breaks = c( ymin, 0.00, ymax))+
             scale_x_continuous(limits = c(0,0.4), breaks = c(0.1,0.3))+
             theme(axis.line =element_line(color = 'black',linewidth=0.3, lineend = 'round') ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.3, lineend = 'round'),
                   axis.ticks.length=unit(.04, "cm"),
                   axis.title =element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(size=5 ,vjust = -7, hjust = 0.15),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0,0.1,0.25,0), 'lines') ) )
  
}

## Calculate the number of pairwise differences in slope per mutation 
# https://strengejacke.github.io/ggeffects/articles/introduction_comparisons.html#hypothesis-testing-for-slopes-of-numeric-predictors

getPropSigDifSlopes = function(data, pCutoff = 0.05)
{
  sub = data #subset(df_mutsCanFitLine, Mut_ID == m)
  sub$Mean_GR = as.numeric(sub$Mean_GR)
  m.interaction <- lm(avgS ~ Mean_GR*env, data =sub)
  
  
  # slopes_withCI = hypothesis_test(m.interaction, c("Mean_GR", "env"), test = NULL) # FOR 71, all but 37sc5 are sig
  slopedifftest = hypothesis_test(m.interaction, c("Mean_GR", "env"))
  slopedifftest$p.value_adj =   p.adjust(slopedifftest$p.value, method = 'BH') # slopedifftest$p.value
  numSigDiffs = sum(slopedifftest$p.value_adj <= pCutoff)
  numTests = length(slopedifftest$p.value_adj <= pCutoff)
  
  out = as.data.frame(slopedifftest)
  sigEnvpairs = paste(subset(out, p.value_adj<=pCutoff)$env, collapse=':')
  alltestedEnvPairs = paste(out$env, collapse=':')
  
  return(data.frame(PropSigDiff = numSigDiffs/numTests, sigEnvpairs = sigEnvpairs, alltestedEnvPairs =alltestedEnvPairs ))
  
}

allPairwiseSlopeDiffs = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(PropSigDiffs_out = map(data, ~getPropSigDifSlopes(.)))%>%
  unnest(PropSigDiffs_out)

allPairwiseSlopeDiffs = allPairwiseSlopeDiffs[,c('Mut_ID', 'PropSigDiff', 'sigEnvpairs','alltestedEnvPairs' )]


#


## Get single slope fit per mutation ##
perMut_1Slope6int_slope = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy)

perMut_1Slope6int_slope = subset(perMut_1Slope6int_slope, term == 'Mean_GR')# only need slopes, not ints
perMut_1Slope6int_slope = perMut_1Slope6int_slope[, c('Mut_ID', 'estimate')]
perMut_1Slope6int_slope$rank_bySlope = rank(perMut_1Slope6int_slope$estimate, ties.method = 'first')# just make the early mutID be 1st rank when are same, since dont want multiple at same place



# Make all main plots
joinedRanks = inner_join(allPairwiseSlopeDiffs,perMut_1Slope6int_slope, by = 'Mut_ID')
joinedRanks$rank_bySlope = factor(joinedRanks$rank_bySlope )
joinedRanks = joinedRanks[(order(joinedRanks$rank_bySlope )),] # put in slope order, so when bin by propSigdiff, will be in slope order
joinedRanks$rowID = 1:nrow(joinedRanks)
allPlots = map(joinedRanks$rowID, ~plotExample_includepairwiseSlopeDiff( .x))
#allPlots_slopeOrder = ggarrange(plotlist = allPlots, nrow = 12, ncol = 8)


## START : Make the insets - half tile plot ###
envs = unique(df_mutsCanFitLine$env)
env1 = numeric(15)
env2 = numeric(15)
counter = 0
for(e1i in 1:(length(envs)-1))
{
  e1 = envs[e1i]
  for(e2i in (e1i+1):length(envs))
  {
    counter = counter + 1
    e2 = envs[e2i]
    env1[counter] = e1
    env2[counter] = e2
  }
}

allEnvPairs_df =  data.frame(env1,env2)
allEnvPairs_df$envPair =  paste0(allEnvPairs_df$env1,'-', allEnvPairs_df$env2)
allEnvPairs_df$altenvPair =  paste0(allEnvPairs_df$env2,'-', allEnvPairs_df$env1)

# flip right 
allEnvPairs_df = allEnvPairs_df[rev(order((allEnvPairs_df$env1))), ]
allEnvPairs_df$env1 = factor(allEnvPairs_df$env1, levels = unique(allEnvPairs_df$env1) )


makeInsetTile = function(rowID)
{
  pairs = joinedRanks$sigEnvpairs[rowID]
  indivSig = str_split(pairs, ":")[[1]]
  allTested = str_split( joinedRanks$alltestedEnvPairs[rowID], ":")[[1]]
  allEnvPairs_df$sigThisMut = allEnvPairs_df$envPair %in% indivSig | allEnvPairs_df$altenvPair %in% indivSig
  allEnvPairs_df$testedThisMut = allEnvPairs_df$envPair %in% allTested | allEnvPairs_df$altenvPair %in% allTested
  allEnvPairs_df$testedThisMut = factor(allEnvPairs_df$testedThisMut, levels = c(T, F))
  
  #  filling ones that are not sig but were tested with grey
  allEnvPairs_df$sigThisMut[allEnvPairs_df$sigThisMut == F & allEnvPairs_df$testedThisMut == T]  = 'realFalse'
  allEnvPairs_df$sigThisMut = factor(allEnvPairs_df$sigThisMut, levels = c(F,T, 'realFalse'))
  allEnvPairs_df$annotate = ''
  allEnvPairs_df$annotate[allEnvPairs_df$sigThisMut == F] = 'X'
  
  plot_tile = ggplot(allEnvPairs_df, aes(x = env1, y = env2, fill = sigThisMut))+
    geom_tile( color = '#948982', linewidth = 0.2)+ # 
    scale_fill_manual(values = c('grey', 'black','white'), drop = F)+
   #scale_color_manual(values = c('black', 'transparent'), drop = F)+
    #geom_text(aes(label=annotate), color = 'black', size = 0.6)+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme_classic()+
    theme_void()+
    theme(axis.text = element_blank(), axis.title = element_blank(), 
          legend.position = 'none' , axis.line = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0.02,0.02,0,0), 'lines') )
  
  return(plot_tile)
}


# use joined ranks so in same order as all Main plots
allInsets = map(1:nrow(joinedRanks), ~makeInsetTile(.x))

## END : Make the insets - half tile plot ###


## put main plots and insets together 
mainPLot_plusInset = function(i)
{
  p1 = allPlots[[i]]
  p2 = allInsets[[i]]
  
  ## custom y axis scale for plot
  yaxisScales$isIn = i >= yaxisScales$mMin & i <= yaxisScales$mMax
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  return(p1 + annotation_custom(
    ggplotGrob(p2), 
    xmin = 0.28, xmax = 0.4, ymin = 0.02, ymax = ymax+0.01 # 0.1
  ))
}

allMutPLots_w_Inset = map(1:94, ~mainPLot_plusInset(.x))
allPlots_slopeOrder_wInset = ggarrange(plotlist = allMutPLots_w_Inset, nrow = 12, ncol = 8)

## to make it so can go around legend 
sub1 = allMutPLots_w_Inset[1:10]
sub2 = allMutPLots_w_Inset[11:82]
sub3 = allMutPLots_w_Inset[83:94]

allPlots_slopeOrder_wInset_1 = ggarrange(plotlist = sub1, nrow = 2, ncol = 5)
allPlots_slopeOrder_wInset_2 = ggarrange(plotlist = sub2, nrow = 9, ncol = 8)
allPlots_slopeOrder_wInset_3 = ggarrange(plotlist = sub3, nrow = 2, ncol = 6)

# ggsave(paste0(figSave,'allPlots_slopeOrder_wInset_1a.pdf' ),allPlots_slopeOrder_wInset_1, width = 16.5/8*5, height = 16.5/12*2, unit = 'cm')
# ggsave(paste0(figSave,'allPlots_slopeOrder_wInset_2a.pdf' ),allPlots_slopeOrder_wInset_2, width = 16.5/8*8, height = 16.5/12*9, unit = 'cm')
# ggsave(paste0(figSave,'allPlots_slopeOrder_wInset_3a.pdf' ),allPlots_slopeOrder_wInset_3, width = 16.5/8*6, height = 16.5/12*2, unit = 'cm')



## ~~ Plot : Barplot % Sig Diff Slopes ####
cBin1 = sum(allPairwiseSlopeDiffs$PropSigDiff == 0)
cBin2 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0 & allPairwiseSlopeDiffs$PropSigDiff <= 0.2)
cBin3 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.2 & allPairwiseSlopeDiffs$PropSigDiff <= 0.4)
cBin4 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.4 & allPairwiseSlopeDiffs$PropSigDiff <= 0.6)
dfPsigDifCol = data.frame(bin = c(1,2,3,4), count = c(cBin1, cBin2,cBin3,cBin4))
dfPsigDifCol_maintext = dfPsigDifCol
allPairwiseSlopeDiffs_maintext = allPairwiseSlopeDiffs
histPsigDiff_realDat = ggplot(dfPsigDifCol, aes(x = as.factor(bin), y = count))+
  geom_col(fill = 'darkgrey',  linewidth = 0.5)+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(0,30,60))+
  ylab('Count')+
  theme_classic()+
  theme(axis.title.y = element_text(color = 'black', size = 10),axis.text = element_text(color = 'black',size = 8), plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.title.x = element_blank(), 
        axis.line =element_line(color = 'black',linewidth=0.3, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.3, lineend = 'round'))
#ggsave(paste0(figSave,'histPsigDiff_MAINTEXT.pdf' ),histPsigDiff_realDat, width = 5.8, height = 2.4, unit = 'cm')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure 4 : Macroscopic Epistasis (DFEs) ####

# ~~ Analysis : Get F0 values ####
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
# ~~ first, need to get slope and intercepts for the model out

Mut_ID = c()
env = c()
slope =c()
intercept =c()
GR = c()
Rsq  = c()
pval = c()
mean_int = c()
resVar = c()
allResiduals = c()
mcounter = 0
for(M in unique(df_mutsCanFitLine$Mut_ID))
{
  mcounter = mcounter + 1
  print(mcounter)
  subMut = subset(df_mutsCanFitLine, Mut_ID == M)
  
  if(dim(subMut)[1] > minNumStrains & length(unique(subMut$env))>1)
  {
    
    
    # ordinary least squares 
    lmout = (lm(avgS ~ Mean_GR + env, data = subMut, na.action = na.exclude)) 
    allCoeff  = lmout$coefficients
    
    
    # fit1 = lmout
    # test = subMut
    # test$predS = fit1$fitted.values
    # # plot with equal slopes, this only really works for mut with all envs measured
    # plot =  ggplot(test,aes(x = Mean_GR))+
    #   geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed', linewidth = 3)+
    #   geom_point(aes(y = avgS, color = env), alpha = 0.5, size  = 3)+
    #   geom_line(aes(y = predS, color = env), linewidth = 3)+
    #   scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'))+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']], color ='#00026E', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env30SC5']], color ='#295DCC', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env30SC7']], color ='#90B6DD', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC3']], color ='#C21700', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC5']], color ='#FE8D26', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC7']], color ='#FBCF96', size = 7 )+
    #   xlim(0,0.4)+
    #   theme_classic() +
    #   xlab('Growth Rate')+
    #   ylab('Selection Coefficient')+
    #   ggtitle(paste0('Mutation ID: ',M))+
    #   theme(axis.line =element_line(color = 'black',linewidth=1.3) ,axis.ticks  =element_line(color = 'black',linewidth=1.3),
    #         axis.title =element_blank(),
    #         axis.text = element_text(color = 'black', size = 15),
    #         plot.title = element_text(size=20),
    #         legend.position = 'none' )
    
    
    
    
    
    for(envChoose in unique(subMut$env))
    {
      
      
      idName = paste0('env', envChoose)
      
      # calc all real ints
      test = allCoeff[-2] # remove GR term
      test[2:length(test)] = test[2:length(test)]+test[1]
      
      if(idName %in% names(allCoeff) ) # if is one of the actual ones
      {
        Mut_ID = c(Mut_ID, M)
        env = c(env, envChoose)
        slope = c(slope, allCoeff[['Mean_GR']] )
        intercept = c(intercept, allCoeff[['(Intercept)']]+allCoeff[[idName]] )
        Rsq  = c(Rsq, summary(lmout)$r.squared)
        pval = c(pval, lmp(lmout))
        mean_int = c(mean_int, mean(test))
        #resVar = c(resVar, var(summary(lmout)$residuals[subMut$env == envChoose], na.rm = T)) # juet get resisduals for this ENV 
        resVar = c(resVar, var(summary(lmout)$residuals, na.rm = T)) #  get resisduals for all ENVs
        
        allResiduals = c(allResiduals, as.numeric(summary(lmout)$residuals))
      }else{ # is the first one, 
        Mut_ID = c(Mut_ID, M)
        env = c(env, envChoose)
        slope = c(slope, allCoeff[['Mean_GR']] )
        intercept = c(intercept, allCoeff[['(Intercept)']])
        Rsq  = c(Rsq, summary(lmout)$r.squared)
        pval = c(pval, lmp(lmout))
        mean_int = c(mean_int, mean(test))
       # resVar = c(resVar,  var(summary(lmout)$residuals[subMut$env == envChoose], na.rm = T)) # juet get resisduals for this ENV 
         resVar = c(resVar, var(summary(lmout)$residuals, na.rm = T)) # juet get resisduals for this ENV 
        
         allResiduals = c(allResiduals, as.numeric(summary(lmout)$residuals))
      }
      
      
      
      
    }
    
    
    
  }else{
    Mut_ID = c(Mut_ID, M)
    env = c(env, envChoose)
    slope = c(slope, NA)
    intercept = c(intercept, NA)
    Rsq  = c(Rsq, NA)
    pval = c(pval, NA)
    mean_int = c(mean_int, NA)
    resVar = c(resVar,  NA)
   
  }
  
}


df_Mutslopes_1slope_6int = data.frame(Mut_ID ,
                                      env ,
                                      slope ,
                                      intercept,
                                      Rsq  ,
                                      pval ,
                                      mean_int ,
                                      resVar)
df_Mutslopes_constantSlopes = df_Mutslopes_1slope_6int


df_F0_1 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), tidy = map(fit, ~ tidy(.x))) %>% # enforced intercept at 0,0
  unnest(tidy) 
df_F0_1 = df_F0_1[, c('env', 'term' ,  'estimate', 'std.error')]
df_F0 = spread(df_F0_1, term, estimate)
names(df_F0) = c('env',  'stdErr','F0')
df_F0$F0 = abs(df_F0$F0)


#write.csv(df_F0, paste0(figSave, 'df_F0.csv'))



# ~~ Analysis : calculate adjusted GR for all strains in all envs ####

df_mutsCanFitLine_wAdjGR = inner_join(df_mutsCanFitLine, df_F0, by = 'env')
df_mutsCanFitLine_wAdjGR$adjGR = df_mutsCanFitLine_wAdjGR$Mean_GR - df_mutsCanFitLine_wAdjGR$F0




# ~~ Plot : All Regressions each env, eta = 0 and otherwise ####
plotAllRegressions = function(environment, eta_to_0 = F, title = T, useAdjGR  = F, addGRrange = T,xlim = c(0,0.8), xbreaks = c(0,0.4,0.8))
{
  
  subDF = subset(df_Mutslopes_1slope_6int, env == environment)
  F0 = subset(df_F0, env == environment)$F0
  xIntVal = F0
  if(eta_to_0 == T & useAdjGR == F) # if using eta to 0, change intercept to that pred by -1*slope*F0
  {
    subDF$intercept = -1*subDF$slope * F0
  }
  
  if(useAdjGR == T & eta_to_0 == F) 
  {
    subDF$intercept =  subDF$intercept + 1*subDF$slope * F0 # shift all ints by mean int pred
    xlim = c(-0.3,0.3)
    xIntVal = 0
  }
  if(eta_to_0 == T & useAdjGR == T)
  {
    subDF$intercept = 0 
    xlim = c(-0.3,0.3)
    xIntVal = 0
  }
  
  testerdata = data.frame(x = c(0.037,F0,0.38), y = 0.1)
  titleUse = paste0(environment)[title] # gives title if T, blank title if F
  
  if(addGRrange )
  {
    
    ## get gr range
    getGRs = range(subset(df_s_ests_wGR, env == environment)$Mean_GR)
    if(useAdjGR == T ) 
    {
      getGRs = getGRs - F0
    }
    
    plotout = ggplot()+
      geom_abline(data = subDF, aes(slope = slope, intercept = intercept, group = Mut_ID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.2,lineend='round')+
      geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
      geom_vline(xintercept = xIntVal, color = 'grey', linewidth = 0.5, linetype = 'solid')+
      annotate("rect", xmin = getGRs[1], xmax = getGRs[2], ymin = -Inf, ymax = -0.4, alpha = 0.2)+
      scale_color_manual( values = c('black', 'black'))+
      scale_x_continuous(limits = xlim, breaks = xbreaks)+
      scale_y_continuous(limits = c(-0.4,0.2), breaks = c(-0.4,-0.2,0,0.2))+
      xlab('Growth Rate (1/hr)')+ ylab('Selection Coefficient')+
      #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
      ggtitle(paste0(titleUse))+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.25, color = 'black',lineend='round'),
            axis.ticks = element_line(linewidth = 0.25, color = 'black',lineend='round'),
            axis.title = element_blank(),
            axis.text =  element_blank(), #text(color = 'black', size = 8),
            legend.text = element_text(color = 'black', size = 8),
            legend.key.size = unit(0.2, 'in'),
            legend.position = 'none',
            plot.margin = unit(c(0.06,0.4,0,0), 'lines')) #text(color = 'black', size = 10,angle = 90)
    
    
  }else{
    plotout = ggplot()+
      geom_abline(data = subDF, aes(slope = slope, intercept = intercept, group = Mut_ID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.2,lineend='round')+
      geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
      geom_vline(xintercept = xIntVal, color = 'grey', linewidth = 0.5, linetype = 'solid')+
      scale_color_manual( values = c('black', 'black'))+
      scale_x_continuous(limits = xlim)+
      scale_y_continuous(limits = c(-0.4,0.2))+
      xlab('Growth Rate (1/hr)')+ ylab('Selection Coefficient')+
      #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
      ggtitle(paste0(titleUse))+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.5, color = 'black',lineend='round'),
            axis.ticks = element_line(linewidth = 0.5, color = 'black',lineend='round'),
            axis.title = element_blank(),
            axis.text =  element_text(color = 'black', size = 8),
            legend.text = element_text(color = 'black', size = 8),
            legend.key.size = unit(0.2, 'in'),
            legend.position = 'none') #text(color = 'black', size = 10,angle = 90)
    
    
  }
 
  return(plotout)
}



SingleEnv_ex = plotAllRegressions( '30SC7', useAdjGR = T, eta_to_0 = T,title = F)
#ggsave(paste0(figSave,'SingleEnv_ex_allReg.pdf' ),SingleEnv_ex, width = 4.9, height = 4.33+0.287+0.358, unit = 'cm')

allEnvs = map(c('30SC3','30SC5','30SC7','37SC3','37SC5','37SC7'),~plotAllRegressions(.x, useAdjGR = F, eta_to_0 = F,title = F) )
allEnvs_a = ggarrange(plotlist = allEnvs, nrow = 3, ncol = 2)


#ggsave(paste0(figSave,'allEnvs_regs_a.pdf' ),allEnvs_a, width = 5.3, height = 7.3, unit = 'cm')


# ~~ Plot : 3 adj GR bins, plot compiled DFEs each env ####
numBins = 3

## Make ADJ GR bins 
binSeq_adjGR = seq(-0.16, 0.16, length.out = numBins+1 )# slightly over range of actual adj GRs, but allows nice #s
adjGR_binDF = data.frame(minBin = binSeq_adjGR[1:(length(binSeq_adjGR)-1)], maxBin =  binSeq_adjGR[2:(length(binSeq_adjGR))])
adjGR_binDF$binID = 1:nrow(adjGR_binDF)
names(adjGR_binDF ) = paste0(names(adjGR_binDF), '_adjGR')

## Make ABS GR bins
binSeq_absGR = seq(0,0.4, length.out = numBins+1 )
absGR_binDF = data.frame(minBin = binSeq_absGR[1:(length(binSeq_absGR)-1)], maxBin =  binSeq_absGR[2:(length(binSeq_absGR))])
absGR_binDF$binID = 1:nrow(absGR_binDF)
names(absGR_binDF ) = paste0(names(absGR_binDF), '_absGR')

## add to df_mutsCanFitLine_wAdjGR
df_mutsCanFitLine_wAdjGR$bin_absGR = NA
df_mutsCanFitLine_wAdjGR$bin_adjGR = NA
for(i in 1:max(absGR_binDF$binID_absGR))
{
  minBin_abs = absGR_binDF$minBin_absGR[i]
  maxBin_abs = absGR_binDF$maxBin_absGR[i]
  
  minBin_adj = adjGR_binDF$minBin_adjGR[i]
  maxBin_adj = adjGR_binDF$maxBin_adjGR[i]
  
  
  df_mutsCanFitLine_wAdjGR$bin_absGR[df_mutsCanFitLine_wAdjGR$Mean_GR >= minBin_abs & df_mutsCanFitLine_wAdjGR$Mean_GR < maxBin_abs] = i
  df_mutsCanFitLine_wAdjGR$bin_adjGR[df_mutsCanFitLine_wAdjGR$adjGR >= minBin_adj & df_mutsCanFitLine_wAdjGR$adjGR < maxBin_adj] = i
  
}


# if put all envs together
plot_compiledDFE_adjGR_allenv = function( binChoose, binwidthChoose = 0.03, adjGR = T)
{
  if(adjGR) # use adjGR
  {
    data = subset(df_mutsCanFitLine_wAdjGR, bin_adjGR == binChoose)
    numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
    envsChoose = subset(numStrainsperenv, numStrains>=3)$env
    data = subset(data, env%in%envsChoose)
    data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
    
  }else{
    data = subset(df_mutsCanFitLine_wAdjGR, bin_absGR == binChoose)
    numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
    envsChoose = subset(numStrainsperenv, numStrains>=3)$env
    data = subset(data, env%in%envsChoose)
    data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
    
  }
  
  
  
  # use stat_binline instead of staat_bin,, this actually just outlines histogram instead of creating shifted bin count
  adjGR_compDFEs = ggplot(data, aes(x = avgS))+
  # geom_histogram(bins = 17, aes(y = after_stat(density)), fill = '#74706e' , color = '#74706e')+
   # stat_binline(bins = 17, aes(color = env,group = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.002),  linewidth = 0.7, alpha = 1 )+
    geom_histogram(binwidth = binwidthChoose, aes(y = after_stat(density)*binwidthChoose), fill = 'grey' , color = 'grey')+
    stat_binline(binwidth = binwidthChoose, aes(color = env,group = env, y = after_stat(density)*binwidthChoose),geom="step", position = position_dodge(width = 0.002),  linewidth = 0.3, alpha = 1 )+
    geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.7,lineend='round')+
    scale_x_continuous(limits = c(-0.2, 0.2), breaks = c(-0.15, 0, 0.15))+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.8), breaks = c(0,0.3,0.6))+
    scale_color_manual(values = c("#00026E", "#295DCC", "#90B6DD", "#C21700" ,"#FE8D26", "#FBCF96"), drop = F)+
    #ggtitle(paste0('Adjusted GR Bin, ', e2))+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.25,lineend='round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.25,lineend='round'), 
          axis.title = element_blank(),
          axis.text.y = element_text(color = 'black', family = 'Helvetica', size = 8),
          axis.text.x = element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(0.2,0.1,0,0), 'lines'))
  
  
  
  return(adjGR_compDFEs)
}

allBinDFEs = ggarrange(plotlist = map(unique(adjGR_binDF$binID_adjGR), ~plot_compiledDFE_adjGR_allenv( .x)), nrow = 1, ncol = length(unique(adjGR_binDF$binID_adjGR)))
plotlist = map(unique(adjGR_binDF$binID_adjGR), ~plot_compiledDFE_adjGR_allenv( .x))
 # ggsave(paste0(figSave,'compiledDFE_bin1.pdf' ),plotlist[[1]], width = 5.62, height = 1.3, unit = 'cm')
 # ggsave(paste0(figSave,'compiledDFE_bin2.pdf' ),plotlist[[2]], width = 5.62, height = 1.3, unit = 'cm')
 # ggsave(paste0(figSave,'compiledDFE_bin3.pdf' ),plotlist[[3]], width = 5.62, height = 1.3, unit = 'cm')




# get dist of GRs for each of the bins
plot_distGRs_eachBin = function( binChoose)
{
  
  
  ## to highlight the area of GRs in the full hist
 minBin =  subset(adjGR_binDF, binID_adjGR ==binChoose)$minBin_adjGR
 maxBin =  subset(adjGR_binDF, binID_adjGR ==binChoose)$maxBin_adjGR
 
  
  data = df_mutsCanFitLine_wAdjGR # subset(df_mutsCanFitLine_wAdjGR, bin_adjGR == binChoose)
  numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
  envsChoose = subset(numStrainsperenv, numStrains>=3)$env
  data = subset(data, env%in%envsChoose)
  data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
  
  data$colorBy = 'outRange'
  data$colorBy[data$bin_adjGR == binChoose] = 'inRange'
  # 
  # adjGR_compDFEs = ggplot(data, aes(x = adjGR))+
  #   #geom_histogram(fill = 'transparent', color = 'black', aes(y = after_stat(density)), bins = 20)+
  # #  geom_rect(xmin=minBin, xmax=maxBin, ymin=0, ymax=6, fill='grey', color="grey", alpha=0.5)+
  #   geom_vline(xintercept = minBin, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #   geom_vline(xintercept = maxBin, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #  stat_binline(aes(y = after_stat(density)),geom="step",bins=10, linewidth = 0.2, alpha = 1)+
  #  # geom_histogram(  aes(y = after_stat(density), fill = colorBy), color = 'black', bins = 10)+
  #   geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #   scale_x_continuous(limits = c(-0.2,0.2), breaks = c(-0.18, 0, 0.18))+
  #   scale_y_continuous(expand = c(0,0), breaks = c(0,3,6), limits = c(0,6))+
  #   #scale_fill_manual(values = c("white", "darkgrey"), drop = F)+
  #   theme_classic()+
  #   theme(axis.line = element_line(color = 'black',linewidth=0.2, lineend = 'round') ,
  #         axis.ticks  =element_line(color = 'black',linewidth=0.2, lineend = 'round'), 
  #         axis.title = element_blank(),
  #         axis.text = element_blank(),
  #         legend.position = 'none',
  #         plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
  #         axis.ticks.length=unit(.05, "cm"))
  # 
  
  numEntries = length(data$Strain)
  binWidthChoose = 0.0355
  adjGR_compDFEs =  ggplot(data, aes(x = adjGR))+
    geom_histogram(binwidth = binWidthChoose, aes(y = after_stat(count)/numEntries, fill = colorBy), color = 'transparent')+
    stat_binline(binwidth = binWidthChoose, aes(y = after_stat(count)/numEntries),geom="step",linewidth = 0.2, alpha = 1, color = 'black')+
    geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5)+
    scale_fill_manual(values = c('#58595B', 'white'))+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25), breaks = c(0,0.2) )+
    scale_x_continuous( limits = c(-0.2,0.2), breaks = c(-0.2,0,0.2) )+
    ylab('Prop.')+
    xlab('Adj. GR')+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.2, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.2, lineend = 'round'), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x =element_text(size = 4, color = 'black', family = 'Helvetica'),
          legend.position = 'none',
          plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
          axis.ticks.length=unit(.05, "cm"))
  
   
  
  
  return(adjGR_compDFEs)
}

grs_bin1 = plot_distGRs_eachBin(1)
grs_bin2 = plot_distGRs_eachBin(2)
grs_bin3 = plot_distGRs_eachBin(3)
# ggsave(paste0(figSave,'grs_bin1.pdf' ),grs_bin1, width = 1.65, height = 0.7, unit = 'cm')
# ggsave(paste0(figSave,'grs_bin2.pdf' ),grs_bin2,  width = 1.65, height = 0.7, unit = 'cm')
# ggsave(paste0(figSave,'grs_bin3.pdf' ),grs_bin3, width = 1.65, height = 0.7, unit = 'cm')



df_F0_noYPD = df_F0
# ~~ Analysis: Predicted DFEs ####
df_F0_2 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), augment = map(fit, ~ augment(.x))) %>% # enforced intercept at 0,0
  unnest(augment) 

# Parameters needed
Mean_slopeDist = mean(df_Mutslopes_1slope_6int$slope) # can also take avg of avg in each env, yields -0.1507 (almost same)
Var_slopeDist = var(df_Mutslopes_1slope_6int$slope)
Skew_slopeDist = skewness(df_Mutslopes_1slope_6int$slope)

lm_resVar_vs_b = summary(lm(resVar ~0+slope, data =df_Mutslopes_1slope_6int )) # enforce 0,0
alpha_resVar_vs_b = -1*lm_resVar_vs_b$coefficients[1]


# for diff type- where slope indep res var
avgResVar = 0.0002291042 
etaVar = var(df_F0_2$.resid) # e.g., sigma^2 pivot


dfe_parampred = function(GR, environment)
{
  
  F_0_e = df_F0[match(environment,df_F0$env ),]$F0 #subset(df_F0, env == environment)$F0
  # this syntax allows for multiple envs inputted 
  # have to either input 1 or numGR entries though
 
  # for NEW slope-independent res var
  # DFEmean = (GR - F_0_e) * Mean_slopeDist 
  # DFEvar = avgResVar + (GR - F_0_e)^2 * Var_slopeDist + etaVar 
  # DFEskew = (Skew_slopeDist *  (GR - F_0_e)^3)/(( (GR - F_0_e)^2 + (avgResVar+etaVar)/Var_slopeDist)^(3/2)) 
  # 
  
  
  # for NEW slope-dependent res var
  DFEmean = (GR - F_0_e) * Mean_slopeDist 
  DFEvar = -alpha_resVar_vs_b*Mean_slopeDist + etaVar + (GR - F_0_e)^2 * Var_slopeDist  
  DFEskew = (Skew_slopeDist *  (GR - F_0_e)^3 + -3*alpha_resVar_vs_b*(GR - F_0_e)/sqrt(Var_slopeDist))/(( (GR - F_0_e)^2 + (-alpha_resVar_vs_b*Mean_slopeDist +etaVar)/Var_slopeDist)^(3/2)) 
  
  
  

  
  return(data.frame(DFEmean, DFEvar, DFEskew))
  
}


# ~~ Analysis :  get DFE stats ####
oldDFEdat_0 = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err, F0) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T), numMuts = sum(!is.na(avgS)))

oldDFEdat_0$predDFEmean = dfe_parampred(as.numeric(oldDFEdat_0$Mean_GR), oldDFEdat_0$env)$DFEmean
oldDFEdat_0$predDFEvar = dfe_parampred(as.numeric(oldDFEdat_0$Mean_GR), oldDFEdat_0$env)$DFEvar
oldDFEdat_0$predDFEskew = dfe_parampred(as.numeric(oldDFEdat_0$Mean_GR), oldDFEdat_0$env)$DFEskew




# ~~ Analysis : Get Error in all params from Bootstrapping ####
numMutsSamp = 70
numIts = 300
set.seed(NULL)

sampleNmuts_DFEstat = function(data, numMutsSamp)
{
  indecies = 1:nrow(data)
  samps = sample(indecies, numMutsSamp, replace = T)
  samp_s = (data[samps,])$avgS
  DFEmean = mean(samp_s)
  DFEvar = var(samp_s)
  DFEskew = skewness(samp_s)

    return(data.frame(DFEmean,DFEvar,DFEskew) )
  
}

bootError = function(data) # data is subset of df_mutsCanFitLine_wAdjGR for 1 env, 1 strain
{
  # For testing, 
  # data = subset(df_mutsCanFitLine_wAdjGR, Strain == 'LK1-A02'& env == '30SC3')
  if(nrow(data) > 5)
  {
    allDFEmean = numeric(numIts)
    allDFEvar = numeric(numIts)
    allDFEskew = numeric(numIts)
    for(i in 1:numIts)
    {
      out = sampleNmuts_DFEstat(data, numMutsSamp)
      allDFEmean[i] = out$DFEmean
      allDFEvar[i] = out$DFEvar
      allDFEskew[i] = out$DFEskew
    }
    
    return(data.frame(DFEmean_SD  = sd(allDFEmean),DFEvar_SD  = sd(allDFEvar),DFEskew_SD  = sd(allDFEskew)   ))
  }else{
    return(data.frame(DFEmean_SD  =NA,DFEvar_SD  = NA,DFEskew_SD  = NA  ))
  }
  
  
}

allDFEbootSD = df_mutsCanFitLine_wAdjGR %>%
  nest(data = -c(Strain,env)) %>%
  mutate(fit = map(data, ~ bootError(.)) )%>%
  unnest(fit) 

allDFEbootSD = allDFEbootSD[, c('Strain', 'env', 'DFEmean_SD', 'DFEvar_SD', 'DFEskew_SD')]
#write.csv(allDFEbootSD,paste0( figSave, 'allDFEbootSD_300Its.csv'))
#allDFEbootSD = read.csv(paste0( figSave, 'allDFEbootSD_300Its.csv'))


oldDFEdat = full_join(oldDFEdat_0,allDFEbootSD, by = c('Strain', 'env') )
oldDFEdat$DFEmean_min = oldDFEdat$DFEmean - oldDFEdat$DFEmean_SD
oldDFEdat$DFEmean_max = oldDFEdat$DFEmean + oldDFEdat$DFEmean_SD
oldDFEdat$DFEvar_min = oldDFEdat$DFEvar - oldDFEdat$DFEvar_SD
oldDFEdat$DFEvar_max = oldDFEdat$DFEvar + oldDFEdat$DFEvar_SD
oldDFEdat$DFEskew_min = oldDFEdat$DFEskew - oldDFEdat$DFEskew_SD
oldDFEdat$DFEskew_max = oldDFEdat$DFEskew + oldDFEdat$DFEskew_SD

oldDFEdat = oldDFEdat[, !(names(oldDFEdat) %in% c('DFEmean_SD', 'DFEvar_SD', 'DFEskew_SD'))]

## ~~ END: Get Error in all params from Bootstrapping


## ~~~ Add Johnson data 
sub_JohnsonDat = subset(JohnsonDat)# ~~ Use ALL strains from their work , !(Strain %in% unique(oldDFEdat$Strain)))


## ~~ Adjust from relative to raw units 
sub_JohnsonDat$DFEmean = sub_JohnsonDat$DFEmean * 0.9729
sub_JohnsonDat$DFEmean_max = sub_JohnsonDat$DFEmean_max * 0.9729
sub_JohnsonDat$DFEmean_min = sub_JohnsonDat$DFEmean_min * 0.9729


sub_JohnsonDat$DFEvar = sub_JohnsonDat$DFEvar * 0.9729^2
sub_JohnsonDat$DFEvar_min = sub_JohnsonDat$DFEvar_min * 0.9729^2
sub_JohnsonDat$DFEvar_max = sub_JohnsonDat$DFEvar_max * 0.9729^2


ypdF0 = 0.52024718
sub_JohnsonDat$F0 = ypdF0
sub_JohnsonDat$numMuts = NA
df_F0 = rbind(df_F0_noYPD, data.frame(env = 'YPD',stdErr = mean(df_F0$stdErr), F0 = ypdF0))
sub_JohnsonDat$adjGR = sub_JohnsonDat$Mean_GR - sub_JohnsonDat$F0


sub_JohnsonDat$predDFEmean = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEmean
sub_JohnsonDat$predDFEvar = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEvar
sub_JohnsonDat$predDFEskew = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEskew


oldDFEdat = rbind(oldDFEdat,sub_JohnsonDat )


# ~~ Plot : All DFE moments ####

DFEMean_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR,color = env))+
  geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ), color = 'white', shape = 16)+
  geom_point(aes(y = DFEmean ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
#  geom_smooth(aes(y = DFEmean ), method = 'lm',formula = 'y~x', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predDFEmean, color = env, group = env), linewidth = 1, lineend='round')+
  xlab('Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous(limits = c(0,0.8), breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )
#ggsave(paste0(figSave, 'DFEMean_plot_rawGR.pdf'), DFEMean_plot_rawGR, width = 4.29, height = 4.35 )


DFEVar_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR ,color = env))+
  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEvar ), color = 'white', shape = 16)+
  geom_point(aes( y = DFEvar), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar, color = env, group = env), linewidth = 1,lineend='round')+
  xlab('Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous( limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c( 0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )


#ggsave(paste0(figSave, 'DFEVar_plot_rawGR.pdf'), DFEVar_plot_rawGR, width = 4.29, height = 4.35 )



DFESkew_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR,color = env))+
  geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEskew ), color = 'white', shape = 16)+
  geom_point(aes( y = DFEskew ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEskew, color = env, group = env), linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous(limits = c(0,0.8), breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c(-3, 0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )
#ggsave(paste0(figSave, 'DFESkew_plot_rawGR.pdf'), DFESkew_plot_rawGR, width = 4.29, height = 4.35 )


## points to add to adjGr plot
adjGR_binDF$midpoint = (adjGR_binDF$minBin_adjGR + adjGR_binDF$maxBin_adjGR)/2
dfextrapoints = data.frame(x = c(adjGR_binDF$midpoint))
predfestatextrapt = dfe_parampred(dfextrapoints$x + df_F0$F0[df_F0$env == '30SC5'], '30SC5' ) # just turn adj gr back to raw for 1 env (will work for any env)
predfestatextrapt$x = dfextrapoints$x

DFEMean_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR))+
  geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max,color = env), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean,color = env), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean,color = env ), color = 'white', shape = 16)+
  geom_point(aes(y = DFEmean,color = env ), alpha = 0.5, shape = 16)+
      geom_smooth(aes(y = DFEmean ), method = 'lm', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEmean),  color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




DFEVar_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR ,color = env))+
  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ), color = 'white', shape = 16)+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), alpha = 0.5, shape = 16)+
    geom_smooth(aes(y = DFEvar ), method = 'lm',formula = 'y~poly(x,2)', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




## get best fit curve for DFE skew
oldDFEdat$sqGR = oldDFEdat$adjGR^2
lmVar = summary(lm(DFEvar~adjGR+sqGR, data =oldDFEdat ))
DFEvar_p0 = lmVar$coefficients[1,1] * (lmVar$coefficients[1,4]<0.05)
DFEvar_p1 =  lmVar$coefficients[2,1] * (lmVar$coefficients[2,4]<0.05)
DFEvar_p2 = lmVar$coefficients[3,1] * (lmVar$coefficients[3,4]<0.05)
DFEvar_p0_se =  lmVar$coefficients[1,2] 
DFEvar_p1_se =  lmVar$coefficients[2,2]
DFEvar_p2_se =  lmVar$coefficients[3,2]



oldDFEdat$newPredVar = DFEvar_p0 + DFEvar_p1*oldDFEdat$adjGR + DFEvar_p2*oldDFEdat$sqGR
oldDFEdat$cuGR = oldDFEdat$adjGR^3
oldDFEdat$cuGR_scaled = oldDFEdat$cuGR / (oldDFEdat$newPredVar + oldDFEdat$sqGR)^(3/2) 
oldDFEdat$sqGR_scaled = oldDFEdat$sqGR / (oldDFEdat$newPredVar + oldDFEdat$sqGR)^(3/2) 
oldDFEdat$GR_scaled = oldDFEdat$adjGR /(oldDFEdat$newPredVar+ oldDFEdat$sqGR)^(3/2) 


lmSkew = summary(lm(DFEskew~GR_scaled+sqGR_scaled+cuGR_scaled, data =oldDFEdat ))
DFE_3_p0 = lmSkew$coefficients[1,1] * (lmSkew$coefficients[1,4]<0.05)
DFE_3_p1 =  lmSkew$coefficients[2,1] * (lmSkew$coefficients[2,4]<0.05)
DFE_3_p2 = lmSkew$coefficients[3,1] * (lmSkew$coefficients[3,4]<0.05)
DFE_3_p3 = lmSkew$coefficients[4,1] * (lmSkew$coefficients[4,4]<0.05)
DFE_3_p0_se =  lmSkew$coefficients[1,2] 
DFE_3_p1_se =  lmSkew$coefficients[2,2]
DFE_3_p2_se =  lmSkew$coefficients[3,2]
DFE_3_p3_se =  lmSkew$coefficients[4,2]


oldDFEdat$predSkew_new = DFE_3_p0 + DFE_3_p1*oldDFEdat$GR_scaled + DFE_3_p2*oldDFEdat$sqGR_scaled + DFE_3_p3*oldDFEdat$cuGR_scaled


DFESkew_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR,color = env))+
  geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ), color = 'white', shape = 16)+
  geom_point(aes( y = DFEskew ), alpha = 0.5, shape = 16)+
    geom_line(aes(y = predSkew_new),linetype = 'dashed',lineend='round', color = 'black', linewidth = 1)+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )



## local rep
x_adj = 0.25--0.2
x_raw = 0.8

width_Mean_GR = 7.1
adjGR_width = width_Mean_GR*(x_adj/x_raw) + 0.2

# 
# ggsave(paste0(figSave,'DFEMean_plot_adjGR.pdf' ),DFEMean_plot_adjGR, width = adjGR_width , height =4.37-0.4, unit = 'cm')
# ggsave(paste0(figSave,'DFEvar_plot_adjGR.pdf' ),DFEVar_plot_adjGR, width = adjGR_width, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(figSave,'DFEskew_plot_adjGR.pdf' ),DFESkew_plot_adjGR, width = adjGR_width, height = 4.37-0.4, unit = 'cm')

# ggsave(paste0(figSave,'DFEMean_plot_rawGR.pdf' ),DFEMean_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(figSave,'DFEvar_plot_rawGR.pdf' ),DFEVar_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(figSave,'DFEskew_plot_rawGR.pdf' ),DFESkew_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
##        S U P P L E M E N T ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
 adjgrowthRate_data = df_mutsCanFitLine_wAdjGR %>%
   group_by(Strain, env) %>%
   summarize(adjGR = mean(adjGR))
growthRate_data = inner_join(growthRate_data,adjgrowthRate_data, by = c('Strain', 'env') )


colorDF = data.frame(env = c("30SC3" ,"30SC5", "30SC7" ,"37SC3", "37SC5" ,"37SC7" ,"YPD"  ), color = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'))

# Figure S1 : Growth Rate Histograms ####
# ~~ Growth Rate Histogram, per Env ####
perEnvGRhist = function(envChoose)
{
  df = subset(growthRate_data, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = Mean_GR,color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 15, aes(color = env,y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(0,0.44), breaks = c(0,0.2,0.4))+
    scale_y_continuous(limits = c(0,15), expand = c(0,0), breaks = c(0,7,14))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allGRhist = map(unique(growthRate_data$env), ~perEnvGRhist(.x))


allGrhist_a = ggarrange(plotlist = allGRhist, nrow = 1, ncol = 6)
#ggsave(paste0(figSave, 'allGrhist_a.pdf'), allGrhist_a, width = 17, height = 5 , unit = 'cm')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S2 : Data Quality Checks ####

# ~~ Plot: Replicate Correlation (Backgrounds and Mutations), Conf exp, and different references ####
  ## In Script processing individual BC s-ests to final s-ests 


# ~~ Plot: Correlation Correction ####
threshNumMuts_forCor = 10 
df_s_ests_wGR_sub = subset(df_mutsCanFitLine_wAdjGR ,  !is.na(avgS))
rawCor = c()
correctedCor = c()
testCor = c()
Mut_ID = c()
env = c()


sigma_sq_bg =   mean(df_s_ests_wGR_sub$Std_err^2)#
sigma_sq_mut =   mean(df_s_ests_wGR_sub$seS^2)#

for(e in unique(df_s_ests_wGR_sub$env))
{
  subEnv = subset(df_s_ests_wGR_sub, env == e)
  print(e)
  for(m in unique(subEnv$Mut_ID))
  {
    subEnv_Mut = subset(subEnv, Mut_ID == m )
    
    if(nrow(subEnv_Mut)>threshNumMuts_forCor)
    {
      
      # ~~~ SK method 
      cov_e_F_s = cov( subEnv_Mut$Mean_GR,subEnv_Mut$avgS)
      # sigma_sq = mean(subEnv_Mut$Std_err^2)#sum(subEnv_Mut$Std_err^2)/length(subEnv_Mut$Std_err)
      # can check that this x2 is relatively close to  meanSq_stdError_s = sum(subEnv_Mut$seS^2)/length(subEnv_Mut$seS)
      
      pre_var_F = var(subEnv_Mut$Mean_GR)
      pre_var_s = var(subEnv_Mut$avgS)
      
      postVar_s = pre_var_F - sigma_sq_bg
      postVar_F = pre_var_s - sigma_sq_bg -  sigma_sq_mut
      
      
      if(postVar_s<0 | postVar_F < 0)#post_corrected_var_dif_obs_exp<0 | post_corrected_var_predFit<0)
      {
        rawCor = c(rawCor,  cor( subEnv_Mut$Mean_GR,subEnv_Mut$avgS))
        correctedCor = c(correctedCor, 'corrected_varBelow0')
        Mut_ID = c(Mut_ID, m)
        env = c(env, e)
      }else{
        
        post_corrected_cor = (cov_e_F_s + sigma_sq_bg)/(sqrt((postVar_s)*(postVar_F)))
        
        
        rawCor = c(rawCor,  cor( subEnv_Mut$Mean_GR,subEnv_Mut$avgS))
        correctedCor = c(correctedCor, post_corrected_cor)
        Mut_ID = c(Mut_ID, m)
        env = c(env, e)
      }
      
      
      
    }else{
      
      rawCor = c(rawCor, 'belowPointThresh')
      correctedCor = c(correctedCor, 'belowPointThresh')
      Mut_ID = c(Mut_ID, m)
      env = c(env, e)
    }
    
    
  }
  
}


df_corCorrection = data.frame(env, Mut_ID, rawCor, correctedCor)

## look at the actual and corrected corrs

sub_df_corCorrection = subset(df_corCorrection,!(rawCor %in% 'belowPointThresh') &  !(correctedCor %in% c('belowPointThresh', 'corrected_varBelow0'))) # remove the character values 
sub_df_corCorrection$rawCor = as.numeric(sub_df_corCorrection$rawCor)
sub_df_corCorrection$correctedCor = as.numeric(sub_df_corCorrection$correctedCor)

CorrCorrelation_Plot = ggplot(sub_df_corCorrection, aes(x = rawCor, y = correctedCor, color = env))+
  geom_point(alpha = 0.5,shape =16)+
  xlab('Raw Correlation')+ylab('Corrected Correlation')+
 # scale_x_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_abline(slope = 1, intercept = 0, color = 'grey', linewidth = 0.5)+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.45, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.45, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 9.48),
        axis.text =  element_text(color = 'black', size = 9.48),
        legend.text = element_text(color = 'black', size = 9.48),
        legend.key.size = unit(0.6, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


#ggsave(paste0(figSave,'CorrCorrelation_Plot.pdf' ),CorrCorrelation_Plot,unit = 'cm', width = 4.53,height = 5.8 )


# ~~ Plot: Lambda_+ vs Lambda_- ####


allSlopes_deltaGR = df_mutsCanFitLine_wAdjGR %>%
  nest(data = -c(env, Mut_ID))%>%
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR , .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) %>%
  dplyr::select(env, Mut_ID, term, estimate, p.value)%>%
  pivot_wider(names_from = term, 
              values_from = c(estimate,p.value))
names(allSlopes_deltaGR) = c('env', 'Mut_ID','intercept', 'slope_delta' , 'pval_intercept', 'pval_slope')
allSlopes_deltaGR$adjP = p.adjust( allSlopes_deltaGR$pval_slope, method = 'BH')
allSlopes_deltaGR$sig = allSlopes_deltaGR$adjP<0.05
#sum(allSlopes_deltaGR$adjP <0.05)/length(allSlopes_deltaGR$env) # 46%
## same as before (250)




allSlopes_lambdaPlus_minus = allExpGR_data %>%
  nest(data = -c(env, Mut_ID))%>%
  mutate(fit = map(data, ~ lm(avgGR ~ Mean_GR , .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) %>%
  dplyr::select(env, Mut_ID, term, estimate, std.error)%>%
  pivot_wider(names_from = term, 
              values_from = c(estimate,std.error))

names(allSlopes_lambdaPlus_minus) = c('env', 'Mut_ID','intercept', 'slope' , 'stdErr_intercept', 'stdErr_slope')

# for values lower than 1
allSlopes_lambdaPlus_minus = subset(allSlopes_lambdaPlus_minus, !is.na(slope))
allSlopes_lambdaPlus_minus$Pnorm = NA
allSlopes_lambdaPlus_minus$Pnorm[allSlopes_lambdaPlus_minus$slope<1]  = pnorm(allSlopes_lambdaPlus_minus$slope[allSlopes_lambdaPlus_minus$slope<1], 1, allSlopes_lambdaPlus_minus$stdErr_slope[allSlopes_lambdaPlus_minus$slope<1])
allSlopes_lambdaPlus_minus$Pnorm[allSlopes_lambdaPlus_minus$slope>1]  = 1- pnorm(allSlopes_lambdaPlus_minus$slope[allSlopes_lambdaPlus_minus$slope>1], 1, allSlopes_lambdaPlus_minus$stdErr_slope[allSlopes_lambdaPlus_minus$slope>1])
allSlopes_lambdaPlus_minus$Pnorm = allSlopes_lambdaPlus_minus$Pnorm*2
# to account for two sided test

allSlopes_lambdaPlus_minus$adjP = p.adjust( allSlopes_lambdaPlus_minus$Pnorm, method = 'BH')
allSlopes_lambdaPlus_minus$sig = allSlopes_lambdaPlus_minus$adjP<0.05
#sum(allSlopes_lambdaPlus_minus$sig, na.rm = T)/sum(!is.na(allSlopes_lambdaPlus_minus$sig)) # 47%



allSlopes_lambdaPlus_minus$slope_minus1 =  allSlopes_lambdaPlus_minus$slope-1



join = inner_join(allSlopes_lambdaPlus_minus, allSlopes_deltaGR, by = c('env', 'Mut_ID') )

join$sigID =as.character(paste0(as.numeric(join$sig.x), as.numeric(join$sig.y)))



join$sigID[as.character(join$sigID) == '00'] = 'neither'
join$sigID[as.character(join$sigID) == '01'] = 'onlyDelta'
join$sigID[as.character(join$sigID) == '10'] = 'onlyGR'
join$sigID[as.character(join$sigID) == '11'] = 'both'

p1_VarSlopes = ggplot(join, aes(x = slope_minus1, y = slope_delta, shape = sigID, color = sigID, fill = sigID))+
  geom_point(alpha = 0.6)+
  scale_shape_manual(values = c(16, 1, 24,25))+
  scale_color_manual(values = c('#978753', '#96533E', '#13EBB5', '#3E6B60'))+
  scale_fill_manual(values = c('#978753', 'transparent', '#13EBB5', '#3E6B60'))+
  xlab('From: Lambda + vs -')+
  ylab('From: deltaLambda')+
  #  ggtitle('Var Slopes')+
  scale_y_continuous(limits = c(-1.2,0.6), breaks = c(-0.6,0,0.6))+
  scale_x_continuous(limits = c(-1.2,0.6), breaks = c(-0.6,0,0.6))+
  geom_abline(slope = 1, intercept = 0, color = 'grey')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10),
        axis.title.y = element_blank(),
        axis.text =  element_text(color = 'black', size = 8),
        legend.text = element_text(color = 'black', size = 10),
        legend.key.size = unit(0.6, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)




## constant slopes 
allSlopes_deltaGR = df_mutsCanFitLine_wAdjGR %>%
  nest(data = -c( Mut_ID))%>%
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR +env, .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) %>%
  dplyr::select( Mut_ID, term, estimate, p.value)%>%
  filter(term == 'Mean_GR') 

names(allSlopes_deltaGR) = c( 'Mut_ID','term', 'slope_delta' , 'pval_delta')
allSlopes_deltaGR$adjP = p.adjust( allSlopes_deltaGR$pval_delta, method = 'BH')
allSlopes_deltaGR$sig = allSlopes_deltaGR$adjP<0.05
#sum(allSlopes_deltaGR$sig)/length(allSlopes_deltaGR$slope_delta) # 89%


allSlopes_lambdaPlus_minus = allExpGR_data %>%
  nest(data = -c( Mut_ID))%>%
  mutate(fit = map(data, ~ lm(avgGR ~ Mean_GR +env, .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) %>%
  dplyr::select( Mut_ID, term, estimate, std.error)%>%
  filter(term == 'Mean_GR') 

names(allSlopes_lambdaPlus_minus) = c( 'Mut_ID','term', 'slope', 'stdErr_slope' )

# for values lower than 1
allSlopes_lambdaPlus_minus = subset(allSlopes_lambdaPlus_minus, !is.na(slope))
allSlopes_lambdaPlus_minus$Pnorm = NA
allSlopes_lambdaPlus_minus$Pnorm[allSlopes_lambdaPlus_minus$slope<1]  = pnorm(allSlopes_lambdaPlus_minus$slope[allSlopes_lambdaPlus_minus$slope<1], 1, allSlopes_lambdaPlus_minus$stdErr_slope[allSlopes_lambdaPlus_minus$slope<1])
allSlopes_lambdaPlus_minus$Pnorm[allSlopes_lambdaPlus_minus$slope>1]  = 1- pnorm(allSlopes_lambdaPlus_minus$slope[allSlopes_lambdaPlus_minus$slope>1], 1, allSlopes_lambdaPlus_minus$stdErr_slope[allSlopes_lambdaPlus_minus$slope>1])
allSlopes_lambdaPlus_minus$Pnorm = allSlopes_lambdaPlus_minus$Pnorm*2
# to account for two sided test

allSlopes_lambdaPlus_minus$adjP = p.adjust( allSlopes_lambdaPlus_minus$Pnorm, method = 'BH')
allSlopes_lambdaPlus_minus$sig = allSlopes_lambdaPlus_minus$adjP<0.05
#sum(allSlopes_lambdaPlus_minus$sig, na.rm = T)/sum(!is.na(allSlopes_lambdaPlus_minus$sig)) # 90%


allSlopes_lambdaPlus_minus$slope_minus1 =  allSlopes_lambdaPlus_minus$slope-1



join_const = inner_join(allSlopes_lambdaPlus_minus, allSlopes_deltaGR, by = c( 'Mut_ID') )


join_const$sigID =as.character(paste0(as.numeric(join_const$sig.x), as.numeric(join_const$sig.y)))



join_const$sigID[as.character(join_const$sigID) == '00'] = 'neither'
join_const$sigID[as.character(join_const$sigID) == '01'] = 'onlyDelta'
join_const$sigID[as.character(join_const$sigID) == '10'] = 'onlyGR'
join_const$sigID[as.character(join_const$sigID) == '11'] = 'both'
join_const$sigID = factor(join_const$sigID, levels = c('both', 'neither', 'onlyDelta', 'onlyGR'))

p2_ConstSlopes =ggplot(join_const, aes(x = slope_minus1, y = slope_delta, shape = sigID, color = sigID, fill = sigID))+
  geom_point(alpha = 0.8)+
  scale_shape_manual(values = c(16, 1, 24,25), drop = F)+
  scale_color_manual(values = c('#978753', '#96533E', '#13EBB5', '#3E6B60'), drop = F)+
  scale_fill_manual(values = c('#978753', 'transparent', '#13EBB5', '#3E6B60'), drop = F)+
  geom_abline(slope = 1, intercept = 0, color = 'grey')+
  xlab('From: Lambda + vs -')+
  ylab('From: deltaLambda')+
  scale_y_continuous(limits = c(-0.6,0.1), breaks = c(-0.6,-0.3,0))+
  scale_x_continuous(limits = c(-0.6,0.1), breaks = c(-0.6,-0.3,0))+
  # ggtitle('Const Slopes')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10),
        axis.title.y = element_blank(),
        axis.text =  element_text(color = 'black', size = 8),
        legend.text = element_text(color = 'black', size = 10),
        legend.key.size = unit(0.6, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


# ggsave(paste0(figSave, 'diffRegTypes_corr_varSlopes.pdf'), p1_VarSlopes, unit = 'cm', width = 3, height = 5  )
# ggsave(paste0(figSave, 'diffRegTypes_corr_constSlopes.pdf'), p2_ConstSlopes, unit = 'cm', width = 3, height = 5  )




# ~~ Plot: Std_err, per Env ####
perEnvStdErrhist = function(envChoose)
{
  # using growthrate data from proces
  df = subset(growthRate_data, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = log10(Std_err),color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 20, aes(color = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(-4.2,-0.0002), breaks = c(-4,-2,0))+
    scale_y_continuous(limits = c(0,3.3), expand = c(0,0), breaks = c(0,3))+
    # scale_x_continuous(limits = c(-0.007,0.06), breaks = c(0,0.05))+
    #scale_y_continuous(limits = c(0,250), expand = c(0,0), breaks = c(0,100,200))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allSTDEhist = map(unique(growthRate_data$env), ~perEnvStdErrhist(.x))


allSTDEhist_a = ggarrange(plotlist = allSTDEhist, nrow = 1, ncol = 6)

#ggsave(paste0(figSave, 'allSTDEhist_a.pdf'), allSTDEhist_a, width = 17, height = 3, unit = 'cm')


# ~~ Plot: Std_err, per Env ####

perEnvMutSEhist = function(envChoose)
{
  df = subset(df_mutsCanFitLine_wAdjGR, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = log10(seS),color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 20, aes(color = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(-4.2,-0.0002), breaks = c(-4,-2,0))+
    scale_y_continuous(limits = c(0,3.3), expand = c(0,0), breaks = c(0,3))+
    #  scale_x_continuous(limits = c(-0.007,0.06), breaks = c(0,0.05))+
    #scale_y_continuous(limits = c(0,250), expand = c(0,0), breaks = c(0,100,200))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allSTDE_mutshist = map(unique(growthRate_data$env), ~perEnvMutSEhist(.x))


allSTDE_mutshist_a = ggarrange(plotlist = allSTDE_mutshist, nrow = 1, ncol = 6)

#ggsave(paste0(figSave, 'allSTDE_mutshist_a.pdf'), allSTDE_mutshist_a, width = 17, height = 3 , unit = 'cm')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S3 : GxG, GxE interactions, %var in sign####

# ~~ Plot: GxG Effect Re-assortment ####

dfUse =  subset(df_mutsCanFitLine_wAdjGR ,  !is.na(avgS))

countMutSigns_GxG = dfUse %>%
  group_by(Mut_ID,env)%>%
  summarize(numStrainsBen_in = sum(mutCall == 'Beneficial' ),
            numStrainsDelt_in = sum(mutCall == 'Deleterious' ), 
            numStrainsNeutral_in = sum(mutCall == 'Neutral' ), 
            numStrainsIn = length(mutCall))

test_GxG = countMutSigns_GxG %>%
  group_by(env) %>%
  summarize(percent_alwaysBen = sum(numStrainsBen_in == numStrainsIn)/length(numStrainsBen_in),
            percent_alwaysDelt = sum(numStrainsDelt_in == numStrainsIn)/length(numStrainsBen_in),
            percent_alwaysNeutral = sum(numStrainsNeutral_in == numStrainsIn)/length(numStrainsBen_in),
            percent_ben_and_delt = sum(numStrainsBen_in != 0 & numStrainsDelt_in != 0 )/length(numStrainsBen_in)
            , nMut = length(numStrainsBen_in))


test_GxG$other = 1 - test_GxG$percent_alwaysBen -  test_GxG$percent_alwaysDelt -   test_GxG$percent_alwaysNeutral - test_GxG$percent_ben_and_delt

test_GxG = test_GxG[, c('env', 'percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt', 'other')]

test_GxG_g = gather(test_GxG, 'test', 'proportion',percent_alwaysBen:other )

test_GxG_g$test = factor(test_GxG_g$test, levels = c('percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral','other', 'percent_ben_and_delt'))

propMutTypes_gxG = ggplot(test_GxG_g, aes(x = env, y = proportion, fill = test))+
  geom_bar(position = 'stack',stat="identity")+
  # scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#1a9b67','#9b1a4e','#6DB3F5','grey', '#e3ac5b'), drop = FALSE)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1.0))+
  xlab('Environment')+
  ylab('Proportion of Mutations')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line =element_line(color = 'black',linewidth =0.5, lineend = 'round') ,axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.text.y = element_text(color = 'black', size = 10) ,
        axis.text.x =  element_blank(),
        legend.position = 'none' )

#ggsave(paste0(figSave ,'propMutTypes_gxG.pdf'), plot  = propMutTypes_gxG, width = 6.35, height =5, unit = 'cm')



# ~~ Plot: GxE Effect Re-assortment ####

countMutSigns_GxE = dfUse %>%
  group_by(Mut_ID,Strain)%>%
  summarize(numEnvsBen_in = sum(mutCall == 'Beneficial' ),
            numEnvsDelt_in = sum(mutCall == 'Deleterious' ), 
            numEnvsNeutral_in = sum(mutCall == 'Neutral' ), 
            numEnvsIn = length(mutCall))


test_GxE = countMutSigns_GxE %>%
  group_by(Strain) %>%
  summarize(percent_alwaysBen = sum(numEnvsBen_in == numEnvsIn)/length(numEnvsBen_in),
            percent_alwaysDelt = sum(numEnvsDelt_in == numEnvsIn)/length(numEnvsBen_in),
            percent_alwaysNeutral = sum(numEnvsNeutral_in == numEnvsIn)/length(numEnvsBen_in),
            percent_ben_and_delt = sum(numEnvsBen_in != 0 & numEnvsDelt_in != 0 )/length(numEnvsBen_in), 
            nmut = length(numEnvsBen_in))



test_GxE$other = 1 - test_GxE$percent_alwaysBen -  test_GxE$percent_alwaysDelt -   test_GxE$percent_alwaysNeutral - test_GxE$percent_ben_and_delt

test_GxE = test_GxE[, c('Strain', 'percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt', 'other')]

test_GxE_g = gather(test_GxE, 'test', 'proportion',percent_alwaysBen:other )

test_GxE_g$test = factor(test_GxE_g$test, levels = c('other','percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt'))


test_GxE_g = test_GxE_g %>%
  group_by(Strain) %>%
  mutate(totProp_BenDelt = sum(proportion[test == 'percent_ben_and_delt']))

test_GxE_g = test_GxE_g[rev(order(test_GxE_g$totProp_BenDelt )), ]

test_GxE_g$test = factor(test_GxE_g$test, levels = c('percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral','other', 'percent_ben_and_delt'))

test_GxE_g$Strain = factor(test_GxE_g$Strain, levels = unique(test_GxE_g$Strain))


propMutTypes_GxE = ggplot(test_GxE_g, aes(x = Strain, y = proportion, fill = test))+
  geom_bar(position = 'stack',stat="identity")+
  # scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#1a9b67','#9b1a4e','#6DB3F5','grey', '#e3ac5b'), drop = FALSE)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0),breaks = c(0.0,0.5,1.0))+
  xlab('Strain')+
  ylab('Proportion of Mutations')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line =element_line(color = 'black',linewidth =0.5, lineend = 'round') ,axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.text.y = element_blank(),
        axis.text.x =  element_blank(),
        legend.position = 'none' )

#ggsave(paste0(figSave ,'propMutTypes_GxE.pdf'), plot  = propMutTypes_GxE, width = 11.5, height = 5, unit = 'cm')



# ~~ Analysis: Partition Variance for Mutation Sign ####
## Unfortunately, the dataset required for this section is too large to upload to github. It is available via email upon request
df_s_ests_wGR_sub = df_mutsCanFitLine_wAdjGR
df_predMutEffect_FDE_realDat = inner_join(df_s_ests_wGR_sub ,df_Mutslopes_1slope_6int , by = c( 'Mut_ID', 'env') )
df_predMutEffect_FDE_realDat$predAvgS_GlobalEpi =  df_predMutEffect_FDE_realDat$intercept  + ( df_predMutEffect_FDE_realDat$Mean_GR * df_predMutEffect_FDE_realDat$slope)
df_predMutEffect_FDE_realDat$probBen_RealDat = 1 - pnorm(0, mean = df_predMutEffect_FDE_realDat$avgS, sd =   df_predMutEffect_FDE_realDat$seS )

df_predMutEffect_FDE_realDat = subset(df_predMutEffect_FDE_realDat, !(Mut_ID %in% neutralMutIDs))
df_predMutEffect_FDE_realDat$signS_obs = sign(df_predMutEffect_FDE_realDat$avgS)
df_predMutEffect_FDE_realDat$expSign = 2*df_predMutEffect_FDE_realDat$probBen_RealDat - 1
# 
# 
# ## ~~ get all bc s-ests 
# 
# df_allSests_allBCs_wSpreds = inner_join(df_allSests_allBCs, df_predMutEffect_FDE_realDat, by = c('Mut_ID', 'Strain', 'env'))
# df_allSests_allBCs_wSpreds$residual = df_allSests_allBCs_wSpreds$s_est - df_allSests_allBCs_wSpreds$predAvgS_GlobalEpi # obs - expected
# 
# df_GE_w_noise_residuals = df_allSests_allBCs_wSpreds %>%
#   group_by(env, Mut_ID) %>%
#   summarize(resVar_allDat = var(residual, na.rm = T))
# 
# mean(df_predMutEffect_FDE_realDat$resVar) # (0.0002)
# mean(df_GE_w_noise_residuals$resVar_allDat) # order of magnitude higher (0.002)
# 
# df_predMutEffect_FDE_realDat = inner_join(df_predMutEffect_FDE_realDat,df_GE_w_noise_residuals, by = c('env', 'Mut_ID') )
# ## ~~~
# 
# 
# 
# ## do per env ##
# noiseVar = c()
# GE_var = c()
# IDE_var = c()
# env = c()
# plots_varPart = list()
# for(e in unique(df_predMutEffect_FDE_realDat$env))
# {
#   sub = subset(df_predMutEffect_FDE_realDat, !(Mut_ID %in% neutralMutIDs) & env == e )
#   sub$signS_obs = sign(sub$avgS)
#   sub$expSign = 2*sub$probBen_RealDat - 1
#   
#   totVar_obsSign = var(sub$signS_obs )
#   Var_to_noise = (1/(nrow(sub)-1)) * sum((sub$signS_obs  - sub$expSign )^2)
#   
#   
#   ## ~~ Now var attributable to noise + idiosynchro 
#   ## first need estimate of residual variance for regressions --> got for each mean s from regressions in each env (1 slopr 6 int)
#   
#   sub$probBen_globalEpi_wResidError = 1 - pnorm(0, mean = sub$predAvgS_GlobalEpi, sd = sqrt(sub$resVar_allDat)) #sub$seS)
#   sub$expSign_ge = 2*sub$probBen_globalEpi_wResidError - 1
#   Var_to_noise_plus_Ide = (1/(nrow(sub)-1)) * sum((sub$signS_obs  - sub$expSign_ge )^2, na.rm = T)
#   
#   
#   # summary(lm(sub$signS_obs ~sub$expSign_ge))
#   
#   Var_globalEpi = totVar_obsSign - Var_to_noise_plus_Ide
#   Var_IDE = Var_to_noise_plus_Ide - Var_to_noise
#   
#   
#   # proportions
#   noiseVar = c(noiseVar , Var_to_noise/totVar_obsSign )
#   GE_var = c(GE_var,Var_globalEpi/totVar_obsSign)
#   IDE_var = c(IDE_var,Var_IDE/totVar_obsSign)
#   env = c(env, e)
#   
#   
#   dfPLot = data.frame(VarExp  = c(Var_to_noise/totVar_obsSign ,Var_globalEpi/totVar_obsSign,  Var_IDE/totVar_obsSign), Type = c('Noise', 'GE', 'IDE'))
#   
#   plots_varPart[[e]] = ggplot(dfPLot, aes(x = "", y = VarExp, fill = Type)) +
#     geom_col() +
#     geom_text(aes(label = paste0(round(VarExp*100, digits = 1), '%')),
#               position = position_stack(vjust = 0.5), size = 3) +
#     scale_fill_manual(values = c('#6c4832', '#bf977f', '#a99d95'))+
#     coord_polar(theta = "y")+
#     theme_void()+
#     theme(legend.position = 'none')
#   
# }
# 
# df_varInSignExp = data.frame(env,GE_var,IDE_var,noiseVar)
# 
# allvarinsignp = ggarrange(plotlist = plots_varPart, nrow = 1, ncol = 6)
# #ggsave(paste0(figSave, 'allvarinsign.pdf'), allvarinsignp, width = 17.5, height = 7/2, unit = 'cm')
# 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####


# Figure S4 : Best Microscopic Epistasis Model, Differences in Distributions of Slopes and Intercepts ####

# ~~ Analysis: Fitting Alternative Models of Global Epistasis ####

## ~~~ Constant slopes (1 intercept per env, same slope all envs)
perMut_1Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model 

perMut_1Slope6int_fits =perMut_1Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_1Slope6int_fits)[2] = paste0(names(perMut_1Slope6int_fits)[2], '_1s6i')


perMut_1Slope6int_Rsq_p = perMut_1Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_1Slope6int_Rsq_p)[2:5] = paste0(names(perMut_1Slope6int_Rsq_p)[2:5], '_1s6i')


## ~~~ 0 slope 1 int per env
perMut_0Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~  env, .)), 
         glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)

perMut_0Slope6int_fits =perMut_0Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_0Slope6int_fits)[2] = paste0(names(perMut_0Slope6int_fits)[2], '_0s6i')


perMut_0Slope6int_Rsq_p = perMut_0Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_0Slope6int_Rsq_p)[2:5] = paste0(names(perMut_0Slope6int_Rsq_p)[2:5], '_0s6i')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Analysis : Best Model Choice ####

## ~~~ Gather all fits into a list 
df_list <- list(perMut_0Slope6int_fits, perMut_1Slope6int_fits, perMut_6Slope6int_fits)
allFits_EpiModels = df_list %>% purrr::reduce(full_join, by='Mut_ID')

likelihood_ratio_test_0s6i_v_1s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_0s6i[[1]], sub$fit_1s6i[[1]])$'Pr(>Chisq)'[2])
}

likelihood_ratio_test_1s6i_v_6s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_1s6i[[1]], sub$fit_6s6i[[1]])$'Pr(>Chisq)'[2])
}

# Apply the function to each row of the tibble
allFits_EpiModels$rowNum = 1:length(allFits_EpiModels$Mut_ID)
resultLRT_Test <- rowwise(allFits_EpiModels) %>%
  mutate(LRT_p_0s6i_v1s6i =  likelihood_ratio_test_0s6i_v_1s6i(allFits_EpiModels,rowNum),
         LRT_p_1s6i_v6s6i =  likelihood_ratio_test_1s6i_v_6s6i(allFits_EpiModels,rowNum))

resultLRT_Test$padj_0s6i_v1s6i = p.adjust(resultLRT_Test$LRT_p_0s6i_v1s6i, method = 'BH')
resultLRT_Test$padj_1s6i_v6s6i = p.adjust(resultLRT_Test$LRT_p_1s6i_v6s6i, method = 'BH')


# Identify the best model for each mutation 
pCutoff = 0.05
resultLRT_Test$bestModel = 'none'
resultLRT_Test$bestModel[resultLRT_Test$padj_0s6i_v1s6i>pCutoff] = '1s6i' # take 0s6i model as special case of 1 slope 6 intercept model (it is just 1, size 0 slope)
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i>pCutoff] = '1s6i'
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i<pCutoff] = '6s6i'







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~ Test if Pivot GRS sig Diff ####
milo_slopes_ints = milo_slopes_ints[, !(names(milo_slopes_ints) %in% 'X')]
milo_slopes_ints$intercept = 0.9729*milo_slopes_ints$intercept # convert w scaling factor have for GR
milo_slopes_ints$env = 'YPD'
subdf_Mutslopes_1slope_6int = df_Mutslopes_1slope_6int[, names(df_Mutslopes_1slope_6int) %in% names(milo_slopes_ints)]
allSlopes = rbind(subdf_Mutslopes_1slope_6int, milo_slopes_ints)

m.interaction <- lm(intercept ~ slope*env, data =allSlopes)
slopedifftest = hypothesis_test(m.interaction, c("slope", "env"))
slopedifftest$p.value_adj =   p.adjust(slopedifftest$p.value, method = 'BH') # slopedifftest$p.value

slopedifftest = data.frame(slopedifftest)
slopedifftest$env1 = substr(slopedifftest$env, 1,5)
slopedifftest$env2 = substr(slopedifftest$env, 7,11)
slopedifftest$sig = slopedifftest$p.value_adj < 0.05

pivotGRdiffs = ggplot(slopedifftest, aes(x = env1, y = env2, fill = sig))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey', 'black'))+
  geom_text(aes(label=round(p.value_adj,2)), color = 'white', size = 2)+
  #ggtitle('Differences in Pivot GR')+
  xlab('Environment 1')+
  ylab('Environment 2')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')

#ggsave(paste0(figSave, 'pivotGRdiffs.pdf'), pivotGRdiffs, unit = 'cm', width = 6.2, height = 6.75 )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~ Analysis: Stat Differences in Slope and Intercept Distributions ####

df_allMutSlopes_6a_6b = perMut_Env_6Slope6int_slopes_ints_spread
envs = unique(df_allMutSlopes_6a_6b$env)

env1 = c()
env2 = c()
slopePval_KS = c()
intPval_KS = c()
slopePval_t = c()
intPval_t = c()
slopePval_F = c()
intPval_F = c()

for(e1i in 1:(length(envs)-1))
{
  e1 = envs[e1i]
  sub_e1 = subset(df_allMutSlopes_6a_6b, env == e1)
  for(e2i in (e1i+1):length(envs))
  {
    e2 = envs[e2i]
    
    sub_e2 = subset(df_allMutSlopes_6a_6b, env == e2)
    
    # join to only test for same set of muts 
    joined = inner_join(sub_e1,sub_e2, by = c('Mut_ID'))
    
    slopeTest_ks = ks.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_ks = ks.test(joined$intercept.x, joined$intercept.y)
    
    slopeTest_t = t.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_t = t.test(joined$intercept.x, joined$intercept.y)
    
    slopeTest_f = var.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_f = var.test(joined$intercept.x, joined$intercept.y)
    
    
    
    env1 = c(env1, e1)
    env2 = c(env2, e2)
    
    slopePval_KS = c(slopePval_KS, slopeTest_ks$p.value )
    intPval_KS = c(intPval_KS, intTest_ks$p.value)
    slopePval_t = c(slopePval_t,slopeTest_t$p.value)
    intPval_t = c(intPval_t,intTest_t$p.value)
    slopePval_F = c(slopePval_F,slopeTest_f$p.value)
    intPval_F = c(intPval_F,intTest_f$p.value)
    
  }
}


dfSlopeInt_ksTEst = data.frame(
  env1 ,
  env2 ,
  slopePval_KS ,
  intPval_KS ,
  slopePval_t ,
  intPval_t ,
  slopePval_F ,
  intPval_F)

dfSlopeInt_ksTEst$intPval_ks_adj = p.adjust(dfSlopeInt_ksTEst$intPval_KS , method = 'BH')
dfSlopeInt_ksTEst$slopePval_ks_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_KS, method = 'BH')

dfSlopeInt_ksTEst$intPval_t_adj = p.adjust(dfSlopeInt_ksTEst$intPval_t , method = 'BH')
dfSlopeInt_ksTEst$slopePval_t_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_t, method = 'BH')

dfSlopeInt_ksTEst$intPval_f_adj = p.adjust(dfSlopeInt_ksTEst$intPval_F , method = 'BH')
dfSlopeInt_ksTEst$slopePval_f_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_F, method = 'BH')


dfSlopeInt_ksTEst$sigInt_KS =  dfSlopeInt_ksTEst$intPval_ks_adj<0.05
dfSlopeInt_ksTEst$sigSlope_KS = dfSlopeInt_ksTEst$slopePval_ks_adj<0.05 

dfSlopeInt_ksTEst$sigInt_t =  dfSlopeInt_ksTEst$intPval_t_adj<0.05
dfSlopeInt_ksTEst$sigSlope_t = dfSlopeInt_ksTEst$slopePval_t_adj<0.05 


dfSlopeInt_ksTEst$sigInt_f =  dfSlopeInt_ksTEst$intPval_f_adj<0.05
dfSlopeInt_ksTEst$sigSlope_f = dfSlopeInt_ksTEst$slopePval_f_adj<0.05 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot : Stat Differences in Slope and Int Distributions ####

p1_ks = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_KS))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_ks_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Slope - KS-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_ks = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_KS))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_ks_adj, digits = 2) ),color = 'white', size = 2)+
  #  ggtitle('Intercept - KS-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p1_t = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_t))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_t_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Slope - t-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_t = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_t))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_t_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Intercept - t-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')

p1_f = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_f))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_f_adj, digits = 2) ),color = 'white', size = 2)+
  #ggtitle('Slope - F-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_f = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_f))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_f_adj, digits = 2) ),color = 'white', size = 2)+
  #ggtitle('Intercept - F-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')

alltests_slope_int_pairwise = ggarrange(plotlist = list(p1_ks, p1_t,p1_f,p2_ks, p2_t,p2_f), nrow = 2, ncol = 3)
# ggsave(paste0(figSave, 'alltests_slope_int_pairwise.pdf'), alltests_slope_int_pairwise, width = 14.8-5.3,height =6, unit = 'cm' )


# ~~~ Plot of residuals vs gr ####
lm_resVar_vs_b = summary(lm(resVar ~0+slope, data =df_Mutslopes_1slope_6int )) # enforce 0,0
alpha_resVar_vs_b = -1*lm_resVar_vs_b$coefficients[1]

residualVariance_vs_slope = ggplot(df_Mutslopes_1slope_6int, aes(x = slope, y = resVar))+
  geom_point(shape = 16, alpha = 0.5)+
  geom_smooth(method = 'lm', formula = 'y~0+x', se = F, color = 'grey', lineend = 'round', linewidth = 0.5)+
  xlab('Slope')+
  ylab('Residual Variance')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(size=10, color = 'black', family = 'Helvetica'),
        axis.text =  element_text(size=8, color = 'black', family = 'Helvetica'),
        legend.position = 'none') 

#ggsave(paste0(figSave, 'residualVariance_vs_slope.pdf'), residualVariance_vs_slope, width = 5.5, height = 6, unit = 'cm')



# ~~ Analysis : get all eta (e.g residuals of slope int corr) ####
df_F0_2 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), augment = map(fit, ~ augment(.x))) %>% # enforced intercept at 0,0
  unnest(augment) 


# ~~ Plot : all Eta Distributions ####
plotEtaDist = function(envChoose)
{
  sub = subset(df_F0_2, env == envChoose)
  plot = ggplot(sub, aes(x = .resid))+
    geom_histogram(bins =25, alpha =1, position = 'identity', aes(y = after_stat(density)), fill = '#595857')+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.5)+
    stat_function(fun = dnorm, args = list(mean = mean(sub$.resid), sd = sd(sub$.resid)), linewidth = 1)+
    scale_y_continuous(limits = c(0,45), breaks = c(0,20,40))+
    scale_x_continuous(limits = c(-0.062,0.062), breaks = c(-0.05,0,0.05))+
    # ggtitle(paste0(envChoose))+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
          axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
          axis.title = element_blank(), 
          axis.text =  element_blank(),
          legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
          legend.key.size = unit(0.2, 'in'),
          legend.position = 'none') #text(color = 'black', size = 10,angle = 90)
  
  
  return(plot)
}

allEtaDist = map(unique(df_F0_2$env), ~plotEtaDist(.x))
allEtaDist_a = ggarrange(plotlist = allEtaDist, nrow = 2, ncol = 3)
#ggsave(paste0(figSave, 'allEtaDist.pdf'), allEtaDist_a, width = 10.7, height = 5, unit = 'cm')

df_F0_2 %>%
  group_by(env) %>%
  summarize(meanEta = mean(.resid), sdEta = sd(.resid))


## ~~ Get Pivot GR from const vs var slopes model ####
# const is just df_F0

df_F0_1_var = perMut_Env_6Slope6int_slopes_ints_spread %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), tidy = map(fit, ~ tidy(.x))) %>% # enforced intercept at 0,0
  unnest(tidy) 
df_F0_1_var = df_F0_1_var[, c('env', 'term' ,  'estimate', 'std.error')]
df_F0_var = spread(df_F0_1_var, term, estimate)
names(df_F0_var) = c('env',  'stdErr','F0')
df_F0_var$F0 = abs(df_F0_var$F0)


joinF0 = inner_join(df_F0_var, df_F0, by = 'env')

constVvarSlope_pivot = ggplot(joinF0, aes(x = F0.x, y = F0.y, color = env))+
  geom_abline(slope =1, intercept = 0 , color = 'grey', linewidth = 0.5, lineend = 'round')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = FALSE)+
  geom_linerange(aes(xmin = F0.x - stdErr.x, xmax = F0.x + stdErr.x))+
  geom_linerange(aes(ymin = F0.y - stdErr.y, ymax = F0.y + stdErr.y))+
  geom_point(size = 1, shape = 16)+
  xlab('Pivot GR - Variable Slopes')+
  ylab('Pivot GR - Invariable Slopes')+
  xlim(0.1,0.35)+
  ylim(0.1,0.35)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave,'constVvarSlope_pivot.pdf' ),constVvarSlope_pivot, width = 5, height = 5.1, unit = 'cm')




# ~~ Plot: Slope Int dist with Milos data ####

### all regressions
getGRs = range(subset(oldDFEdat, env == 'YPD')$Mean_GR)

allDat_slopes_ints = rbind(perMut_Env_6Slope6int_slopes_ints_spread,milo_slopes_ints)

allDat_slopes_ints$env2 = 'Pooled'
allDat_slopes_ints$env2[allDat_slopes_ints$env == 'YPD'] = 'YPD'


hist_mutSlopes_wYPD = ggplot(data = allDat_slopes_ints, aes(x = slope,color = env2)) + 
  stat_binline(binwidth = 0.12, aes(color = env2,y = after_stat(density)*0.12),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('black', '#6a9c67'), drop = FALSE)+
  scale_fill_manual(values = c('black', '#6a9c67'), drop = FALSE)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.5, lineend = 'round')+
  theme_classic()+
  xlab('Slope')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.position = 'none')
#ggsave(paste0(figSave,'hist_mutSlopes_wYPD.pdf' ),hist_mutSlopes_wYPD,width = 6, height = 5, unit = 'cm')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot:  Slope Int Corr w YPD ####
slopeIntCorr_wYPD = ggplot(allDat_slopes_ints, aes(x = slope, y = intercept, color = env, group = env))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_point( alpha = 0.5, size = 2, shape = 16)+
  geom_smooth(method = 'lm',formula = 'y~0+x', se = F, lineend = 'round')+
  xlab('Slope') +  ylab('Intercept')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = FALSE)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave,'slopeIntCorr_wYPD.pdf' ),slopeIntCorr_wYPD, width = 6, height = 5, unit = 'cm')




# ~~ Plot : F0 vs Mean GR ####

subgr = oldDFEdat %>%
  group_by(Strain, env) %>%
  summarize(Mean_GR = mean(Mean_GR))

perEnvMeanGR = subgr %>%
  group_by(env) %>%
  summarize(Mean_GR_thisEnv = mean(Mean_GR), seGR = sd(Mean_GR)/sqrt(length(Mean_GR)))


df_F0_wMeanGR = inner_join(df_F0,perEnvMeanGR, by = 'env' )
f0_vs_meanGR = ggplot(df_F0_wMeanGR, aes(x =Mean_GR_thisEnv, y = F0, color = env))+
  geom_abline(slope =1, intercept = 0 , color = 'grey', linewidth = 0.5, lineend = 'round')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = FALSE)+
  geom_point(size = 1, shape = 16)+
  geom_linerange(aes(ymin = F0 -stdErr, ymax = F0 +stdErr ), lineend = 'round')+
  geom_linerange(aes(xmin = Mean_GR_thisEnv -seGR, xmax = Mean_GR_thisEnv +seGR ), lineend = 'round')+
  xlab('Average GR')+
  ylab('Pivot GR')+
  xlim(0.1,0.7)+
  ylim(0.1,0.7)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave,'f0_vs_meanGR.pdf' ), f0_vs_meanGR, width = 6, height = 5, unit = 'cm')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####


# Figure S5 :Invariant vs Variant Slopes ####

# ~~ Plot : Best Model ####
barPLot_bestModel = ggplot( resultLRT_Test, aes(x = bestModel))+
  geom_bar(aes(y= after_stat(count)/sum(after_stat(count))*100),fill = 'darkgrey', color = 'darkgrey',  linewidth = 1)+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(breaks = c(0,30,60), expand = c(0,0))+
  theme_classic() +
  xlab('Test')+
  ylab('% Best Model')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'))

#ggsave(paste0(figSave,'barPLot_bestModel.pdf' ),barPLot_bestModel, width = 8.4, height = 5.1-1.5, unit = 'cm')





# ~~ Plot : Violin of Adjusted Rsq  ####
# ~~~ Gather all Rsq and Pvals into a list 
df_list <- list(perMut_0Slope6int_Rsq_p, perMut_1Slope6int_Rsq_p, perMut_6Slope6int_Rsq_p)
allRsq_Pval_EpiModels = df_list %>% purrr::reduce(full_join, by='Mut_ID')

justRsq_EpiModels = allRsq_Pval_EpiModels[,c('Mut_ID', "adj.r.squared_1s6i", "adj.r.squared_6s6i")]
names(justRsq_EpiModels) = c('MutID', '1s6i', '6s6i')
justRsq_EpiModels_g = gather(justRsq_EpiModels, key = 'Model', value = 'Rsq', '1s6i':'6s6i')

deltaRsqPlot = ggplot(data =  justRsq_EpiModels_g, aes(x = Model, y =Rsq)) +
  geom_point( position = 'jitter', alpha = 0.1, shape = 16)+
  #geom_boxplot( colour = 'black', alpha = 0.4, linewidth = 0.4, fill = 'grey') +
  #geom_violin(fill = 'transparent', linewidth = 1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.7,
               colour = "black", lineend = 'round')+
  #scale_fill_manual(values = c('palegreen4','#debe92','#92b2de' ))+
  theme_classic()+
  ylab('Adjusted Rsq')+
  xlab('Model')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8) ,
        legend.position = 'none')

#ggsave(paste0(figSave,'deltaRsqPlot.pdf' ),deltaRsqPlot, width = 4.2, height = 5.1, unit = 'cm')



# ~~ Plot: no missing measures ####


# ~~~ Rsq without missing measures ####

## ~~ Get data with no missing measures per strains 
## e.g., only keep Mut-ID : Strain pairs if that same pair is seen in all envs 
df_s_ests_wGR_countNumEnvIn =  df_s_ests_wGR %>%
  group_by(Mut_ID, Strain) %>%
  mutate(nEnv = length(unique(env)))


sub_mutsInallEnvs = subset(df_s_ests_wGR_countNumEnvIn, nEnv == 6) # this pair in all envs 

sub_mutsInallEnvs =sub_mutsInallEnvs %>%
  group_by(Mut_ID, env) %>%
  mutate(nStrainsforReg = length(Strain)) # down to 90 mutations


minNumStrainsForReg = 5


subTestRegs = subset(sub_mutsInallEnvs, nStrainsforReg>=minNumStrainsForReg)
# nrow(subTestRegs) = 12606, nrow(df_s_ests_wGR ) = 18551
# retaining ~68% of rows


getRsq_1SlopeReg = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)%>%
  dplyr::select(Mut_ID, adj.r.squared, p.value)
getRsq_1SlopeReg$type = 'Constant Slope'


getRsq_1PerEnvSlopeReg = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR * env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)%>%
  dplyr::select(Mut_ID, adj.r.squared, p.value)
getRsq_1PerEnvSlopeReg$type = 'Variable Slope'

joinedRsq_microEpi_noMissing = rbind(getRsq_1SlopeReg,getRsq_1PerEnvSlopeReg )

RsqDist_noMisiingMuts = ggplot(data =  joinedRsq_microEpi_noMissing, aes(x = type, y =adj.r.squared )) +
  geom_point( position = 'jitter', alpha = 0.1, shape = 16)+
  #geom_boxplot( colour = 'black', alpha = 0.4, linewidth = 0.4, fill = 'grey') +
  #geom_violin(fill = 'transparent', linewidth = 1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.7,
               colour = "black",  lineend = 'round')+
  #scale_fill_manual(values = c('palegreen4','#debe92','#92b2de' ))+
  theme_classic()+
  ylab('Adjusted Rsq')+
  xlab('Model')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8) ,
        legend.position = 'none')
#ggsave(paste0(figSave, 'RsqDist_noMisiingMuts.pdf'),RsqDist_noMisiingMuts, unit = 'cm', width = 4.2, height = 5.1-1.5)



# ~~~ Analysis, model fits,  no missing ####

## ~~~ Var slopes
perMut_6Slope6int_Rsq_p_noMissing = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR * env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model 

perMut_6Slope6int_fits_noMissing =perMut_6Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'fit' )]
names(perMut_6Slope6int_fits_noMissing)[2] = paste0(names(perMut_6Slope6int_fits_noMissing)[2], '_6s6i')


perMut_6Slope6int_Rsq_p_noMissing = perMut_6Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_6Slope6int_Rsq_p_noMissing)[2:5] = paste0(names(perMut_6Slope6int_Rsq_p_noMissing)[2:5], '_6s6i')

## ~~~ Constant slopes (1 intercept per env, same slope all envs)
perMut_1Slope6int_Rsq_p_noMissing = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model 

perMut_1Slope6int_fits_noMissing =perMut_1Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'fit' )]
names(perMut_1Slope6int_fits_noMissing)[2] = paste0(names(perMut_1Slope6int_fits_noMissing)[2], '_1s6i')


perMut_1Slope6int_Rsq_p_noMissing = perMut_1Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_1Slope6int_Rsq_p_noMissing)[2:5] = paste0(names(perMut_1Slope6int_Rsq_p_noMissing)[2:5], '_1s6i')


## ~~~ 0 slope 1 int per env
perMut_0Slope6int_Rsq_p_noMissing = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~  env, .)), 
         glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)

perMut_0Slope6int_fits_noMissing =perMut_0Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'fit' )]
names(perMut_0Slope6int_fits_noMissing)[2] = paste0(names(perMut_0Slope6int_fits_noMissing)[2], '_0s6i')


perMut_0Slope6int_Rsq_p_noMissing = perMut_0Slope6int_Rsq_p_noMissing[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_0Slope6int_Rsq_p_noMissing)[2:5] = paste0(names(perMut_0Slope6int_Rsq_p_noMissing)[2:5], '_0s6i')



## ~~~ Gather all fits into a list 
df_list <- list(perMut_0Slope6int_fits_noMissing, perMut_1Slope6int_fits_noMissing, perMut_6Slope6int_fits_noMissing)
allFits_EpiModels = df_list %>% purrr::reduce(full_join, by='Mut_ID')

likelihood_ratio_test_0s6i_v_1s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_0s6i[[1]], sub$fit_1s6i[[1]])$'Pr(>Chisq)'[2])
}

likelihood_ratio_test_1s6i_v_6s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_1s6i[[1]], sub$fit_6s6i[[1]])$'Pr(>Chisq)'[2])
}

# Apply the function to each row of the tibble
allFits_EpiModels$rowNum = 1:length(allFits_EpiModels$Mut_ID)
resultLRT_Test <- rowwise(allFits_EpiModels) %>%
  mutate(LRT_p_0s6i_v1s6i =  likelihood_ratio_test_0s6i_v_1s6i(allFits_EpiModels,rowNum),
         LRT_p_1s6i_v6s6i =  likelihood_ratio_test_1s6i_v_6s6i(allFits_EpiModels,rowNum))

resultLRT_Test$padj_0s6i_v1s6i = p.adjust(resultLRT_Test$LRT_p_0s6i_v1s6i, method = 'BH')
resultLRT_Test$padj_1s6i_v6s6i = p.adjust(resultLRT_Test$LRT_p_1s6i_v6s6i, method = 'BH')

resultLRT_Test = subset(resultLRT_Test, !is.na(resultLRT_Test$padj_0s6i_v1s6i) & !is.na(resultLRT_Test$padj_1s6i_v6s6i))

# Identify the best model for each mutation 
pCutoff = 0.05
resultLRT_Test$bestModel = 'none'
resultLRT_Test$bestModel[resultLRT_Test$padj_0s6i_v1s6i>pCutoff] = '1s6i' # take 0s6i model as special case of 1 slope 6 intercept model (it is just 1, size 0 slope)
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i>pCutoff] = '1s6i'
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i<pCutoff] = '6s6i'




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot : Best Model no missing measure ####
barPLot_bestModel_noMissingMuts = ggplot( resultLRT_Test, aes(x = bestModel))+
  geom_bar(aes(y= after_stat(count)/sum(after_stat(count))*100),fill = 'darkgrey', color = 'darkgrey',  linewidth = 1)+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(breaks = c(0,30,60), expand = c(0,0))+
  theme_classic() +
  xlab('Test')+
  ylab('% Best Model')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'))

#ggsave(paste0(figSave,'barPLot_bestModel_noMissingMuts.pdf' ),barPLot_bestModel_noMissingMuts, width = 8.4, height = 5.1-1.5, unit = 'cm')

sum(resultLRT_Test$bestModel == '1s6i')
sum(resultLRT_Test$bestModel == '6s6i')



# ~~~~~~ prop sig diff slopes no missing ####
getPropSigDifSlopes = function(data, pCutoff = 0.05)
{
  sub = data #subset(df_mutsCanFitLine, Mut_ID == m)
  sub$Mean_GR = as.numeric(sub$Mean_GR)
  m.interaction <- lm(avgS ~ Mean_GR*env, data =sub)
  
  
  # slopes_withCI = hypothesis_test(m.interaction, c("Mean_GR", "env"), test = NULL) # FOR 71, all but 37sc5 are sig
  slopedifftest = hypothesis_test(m.interaction, c("Mean_GR", "env"))
  slopedifftest$p.value_adj =   p.adjust(slopedifftest$p.value, method = 'BH') # slopedifftest$p.value
  numSigDiffs = sum(slopedifftest$p.value_adj <= pCutoff)
  numTests = length(slopedifftest$p.value_adj <= pCutoff)
  
  out = as.data.frame(slopedifftest)
  sigEnvpairs = paste(subset(out, p.value_adj<=pCutoff)$env, collapse=':')
  alltestedEnvPairs = paste(out$env, collapse=':')
  
  return(data.frame(PropSigDiff = numSigDiffs/numTests, sigEnvpairs = sigEnvpairs, alltestedEnvPairs =alltestedEnvPairs ))
  
}

allPairwiseSlopeDiffs = subTestRegs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(PropSigDiffs_out = map(data, ~getPropSigDifSlopes(.)))%>%
  unnest(PropSigDiffs_out)

allPairwiseSlopeDiffs = allPairwiseSlopeDiffs[,c('Mut_ID', 'PropSigDiff', 'sigEnvpairs','alltestedEnvPairs' )]

cBin1 = sum(allPairwiseSlopeDiffs$PropSigDiff == 0)
cBin2 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0 & allPairwiseSlopeDiffs$PropSigDiff <= 0.2)
cBin3 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.2 & allPairwiseSlopeDiffs$PropSigDiff <= 0.4)
cBin4 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.4 & allPairwiseSlopeDiffs$PropSigDiff <= 0.6)
dfPsigDifCol = data.frame(bin = c(1,2,3,4), count = c(cBin1, cBin2,cBin3,cBin4))

histPsigDiff = ggplot(dfPsigDifCol, aes(x = as.factor(bin), y = count))+
  geom_col(fill = 'black', color = '#948982', linewidth = 0.5)+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(0,50), breaks = c(0,25,50))+
  ylab('Count')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),  
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'))
#ggsave(paste0(figSave, 'histPsigDiff_noMissingMuts.pdf'), histPsigDiff, width = 5-0.33, height = 5.1-1.5, unit = 'cm')





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####


# Figure S6 : All Mutation Microscopic Epistasis (w Points) ####


### Plot: all model best fits  with POINTS ###
yaxisScales = data.frame(mMin = c(1,57,81),mMax =  c(56,80,94), scaleMin = c(-0.12,-0.09,-0.05), scaleMax = c(0.12,0.09, 0.05))
# have y-axis scale vary by rows 
plotExample_wPoints = function(rowID) 
{
  m = as.numeric(joinedRanks[rowID,'Mut_ID'])
  
  ## custom y axis scale for plot
  yaxisScales$isIn = rowID >= yaxisScales$mMin & rowID <= yaxisScales$mMax
  ymin = yaxisScales$scaleMin[yaxisScales$isIn]
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  fitted = subset(perMut_6Slope6int_fits, Mut_ID == m)$fit_6s6i[[1]] # get the variable slope model regression
  data = subset(df_mutsCanFitLine, Mut_ID == m)
  data$predS = fitted$fitted.values # note this only works because data here is subsetted from eaxacr same DF in exact same way as when do the LM
  data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  mutName = subset(Mut_full_info_Use, Mut_ID ==m)$Gene.Use
  
  
  #  colorScale_pSigDiff = rev(c('#f8db73', '#f7bae1', '#d42bd6', '#680185'))
  
  #  colorScale_pSigDiff = rev(c('#f8db73', '#d42bd6', '#680185', 'black'))

  return(  ggplot(data,aes(x = Mean_GR))+
             geom_hline(yintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5/2)+
             geom_point(aes(y = avgS, color = env), alpha = 0.2, size  = 0.4, shape = 16)+
             geom_line(aes(y = predS, color = env),   linewidth = 0.5)+
             scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = F)+
             theme_void() +
             xlab('Growth Rate (1/hr)')+
             ylab('Selection Coefficient (1/hr)')+
             # ggtitle(paste0(mutName))+
             annotate(geom="text", x=-Inf, y=-Inf, label=paste0(mutName),
                      color='black', size = 1.945, vjust = -0.6, hjust = -0.1)+
             #scale_y_continuous( breaks = c(-0.12,-0.09, -0.06,-0.03, -0.01, 0.00, 0.01,0.03, 0.06, 0.09, 0.12))+
             scale_y_continuous(limits = c(ymin,ymax), breaks = c( ymin, 0.00, ymax))+
             scale_x_continuous(limits = c(0,0.4), breaks = c(0.1,0.3))+
             theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.5),
                   axis.ticks.length=unit(.04, "cm"),
                   axis.title =element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(size=5 ,vjust = -7, hjust = 0.15),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0,0.1,0.25,0), 'lines') ) )
  
}

allPlots = map(joinedRanks$rowID, ~plotExample_wPoints( .x)) 
## muts 1-56, use -0.09:0.09
## muts 57:80 use -0.05:0.05
## muts 81:94 use - 0.03 : 0.03
# arrange only by slope
allPlots_slopeOrder_wPoints = ggarrange(plotlist = allPlots, nrow = 12, ncol = 8)

#ggsave(paste0(figSave,'allPlots_slopeOrder_wPoints.pdf' ),allPlots_slopeOrder_wPoints, width = 16.5, height = 16.5, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S7 : Growth Rate Re-assortment ####

# ~~ Plot : Correlation of GRs for all Environment Pairs ####
e1Vec = c()
e2Vec = c()
grCorr = c()

allE = unique(growthRate_data$env)
for(e1i in 1:(length(allE)-1))
{
  e1 = allE[e1i]
  sub1 = subset(growthRate_data, env == e1)
  for(e2i in (e1i+1):length(allE))
  {
    e2 = allE[e2i]
    sub2 = subset(growthRate_data, env == e2)
    joined = inner_join(sub1, sub2, by = c('Strain'))
    
    e1Vec = c(e1Vec, e1)
    e2Vec = c(e2Vec, e2)
    grCorr = c(grCorr, cor(joined$Mean_GR.x, joined$Mean_GR.y))
    
  }
}
dfGrCorr = data.frame(e1Vec, e2Vec, grCorr)

allEnvPairs = expand.grid(allE,allE)
names(allEnvPairs) = c('e1Vec', 'e2Vec')

dfGrCorr_wNA = full_join(dfGrCorr, allEnvPairs, by = c('e1Vec', 'e2Vec'))
dfGrCorr_wNA = subset(dfGrCorr_wNA, e1Vec != e2Vec & e1Vec != '37SC7' & e2Vec != '30SC3')

grCorr_tile = ggplot(dfGrCorr_wNA, aes(x = e1Vec, y = e2Vec, fill = grCorr))+
  geom_tile(color = 'white', linewidth = 0.5)+
  geom_text(aes(label=round(grCorr, digits = 2)), color = 'black', size = 3)+
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 0, na.value = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'right')

#ggsave(paste0(figSave ,'grCorr_tile.pdf'), plot  = grCorr_tile, width = 10, height = 7, unit = 'cm')







# ~~ Plot : Growth Rate Re-assortment ####
##~~ add median GR

sumGR = growthRate_data %>%
  group_by(env) %>%
  mutate(medianGR = median(Mean_GR), medianAdjGR = median(adjGR))

sumGR = sumGR[order(sumGR$medianGR), ]
sumGR$env = factor(sumGR$env, levels = unique(sumGR$env))


## stats for changing side
sumGR$topSErange = sumGR$Mean_GR + sumGR$Std_err
sumGR$bottomSErange = sumGR$Mean_GR - sumGR$Std_err

sumGR$belowMedian = sumGR$topSErange < sumGR$medianGR
sumGR$aboveMedian = sumGR$bottomSErange >= sumGR$medianGR

sumGR = sumGR %>%
  group_by(Strain) %>%
  mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)

sumGR$relGR = sumGR$Mean_GR - sumGR$medianGR
sumGR$env = factor(sumGR$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))

relGR_reassort = ggplot(sumGR, aes(x = env))+
  geom_line(aes( y = relGR,group = Strain, color = changeSideofMedian),alpha = 0.5)+
  geom_point(aes( y = relGR, group = Strain, color = changeSideofMedian),alpha = 0.5, shape = 16)+
  scale_color_manual(values = c('black', '#993254'))+#beb1b5
  geom_hline(color = 'grey', linewidth = 0.5, yintercept = 0 )+
  xlab("Environment")+
  ylab('Background Relative Growth Rate (1/h)')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')

#ggsave(paste0(figSave ,'relGR_reassort.pdf'), plot  = relGR_reassort, width = 8, height = 7, unit = 'cm')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S8 : Mutation Effect Re-assortment ####

## first, determine order so can have ones with different shared y axes 
minrange = 0
maxrange = 0
whichAxis = c()
Schoose_all = c()
count = 0
for(Schoose in unique(df_s_ests_wGR$Strain))
{
  count = count + 1
  sub = subset(df_s_ests_wGR, Strain == Schoose)
  
  sub = sub %>%
    group_by(env) %>%
    mutate(medianS = median(avgS))
  
  sub = sub[order(sub$medianS), ]
  sub$env = factor(sub$env, levels = unique(sub$env))
  
  
  ## stats for changing side
  sub$topSErange = sub$avgS + sub$seS
  sub$bottomSErange = sub$Mean_GR - sub$medianS
  
  
  sub$belowMedian = sub$topSErange < sub$medianS
  sub$aboveMedian = sub$bottomSErange >sub$medianS
  
  
  sub = sub %>%
    group_by(Mut_ID) %>%
    mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)
  
  sub$relS = sub$avgS - sub$medianS
  
  minrels = min(sub$relS)
  maxrels = max(sub$relS)
  
  if(minrels > -0.12 & maxrels < 0.12)
  {
    whichAxis = c(whichAxis, '0.06')
  }else{
    whichAxis = c(whichAxis, '0.12')
  }
  
  
  Schoose_all = c(Schoose_all,Schoose )
  
  
}

df_order = data.frame(Schoose_all, whichAxis)
df_order = df_order[order(df_order$whichAxis),]
df_order$Schoose_all = factor(df_order$Schoose_all, levels = unique(df_order$Schoose_all))

## now, make plots
minrange = 0
maxrange = 0
plist_strainMutrearrange = list()
plist_strainMutrearrange_relS = list()
allPchangeSide = c()
count = 0
for(Schoose in unique(df_order$Schoose))
{
  count = count + 1
  sub = subset(df_s_ests_wGR, Strain == Schoose)
  
  sub = sub %>%
    group_by(env) %>%
    mutate(medianS = median(avgS))
  
  sub = sub[order(sub$medianS), ]
  sub$env = factor(sub$env, levels = unique(sub$env))
  
  
  ## stats for changing side
  sub$topSErange = sub$medianS + sub$seS
  sub$bottomSErange = sub$medianS - sub$seS
  
  
  sub$belowMedian = sub$avgS < sub$bottomSErange
  sub$aboveMedian = sub$avgS > sub$topSErange
  
  
  sub = sub %>%
    group_by(Mut_ID) %>%
    mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)
  
  propChangeSide =  sum(sub$changeSideofMedian)/length(sub$changeSideofMedian)
  allPchangeSide = c(allPchangeSide,propChangeSide)
  

  
  sub$relS = sub$avgS - sub$medianS
  sub$env = factor(sub$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  
  minrange = min(minrange, min(sub$relS))
  maxrange = max(maxrange, max(sub$relS))


  plist_strainMutrearrange_relS[[count]] = ggplot(sub, aes(x = env))+
    geom_line(aes( y = relS,group = Mut_ID, color = changeSideofMedian),alpha = 0.5, linewidth = 0.2, lineend = 'round')+
    geom_point(aes( y = relS, group = Mut_ID), color = 'white' , alpha = 0.5, size = 0.35, shape = 16)+ # add underlay for points 
    geom_point(aes( y = relS, group = Mut_ID, color = changeSideofMedian, fill =changeSideofMedian ),alpha = 0.5, size = 0.4, shape = 16)+
    scale_x_discrete(drop = F)+
    scale_y_continuous(limits = c(-0.18, 0.15), breaks = c(-0.12, 0 , 0.12))+
    #scale_y_continuous( breaks = axisUse_breaks)+
    scale_color_manual(values = c('black', '#993254'))+
    geom_hline(yintercept = 0 , color = 'grey', linewidth = 0.5)+
    annotate(geom="text", x=-Inf, y=-Inf, label=paste0(round(propChangeSide, 2)*100, '%'),
             color='#993254', size = 1.945, vjust = -0.6, hjust = -0.1)+
    annotate(geom="text", x=-Inf, y=-Inf, label=paste0(', ', round(1-propChangeSide, 2)*100, '%'),
             color='black', size = 1.945, vjust = -0.6, hjust = -0.9)+
    xlab("Environment")+
    ylab('Relative Fitness Effect (1/h)')+
    ggtitle(paste0(Schoose))+#,subtitle =  paste0(round(propChangeSide, digits = 2)))+
    theme_classic()+
    theme(axis.line =element_line(color = 'black',linewidth=0.3, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.3, lineend = 'round'), 
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(color = 'black',vjust = -4, hjust = 0.1, size = 9, family = 'Helvetica'),
          #plot.subtitle = element_text(color = 'black', size = 7, family = 'Helvetica'),
          legend.position = 'none',
          plot.margin = unit(c(0,0.1,0.25,0), 'lines') )
  
  
  
  
  
  
  
}

allStrains_Mutreassort = ggarrange(plotlist = plist_strainMutrearrange_relS, nrow = 7, ncol = 6)
#ggsave(paste0(figSave ,'allStrains_Mutreassort.pdf'), plot  = allStrains_Mutreassort, width = 16.5, height = 15, unit = 'cm')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####



# Figure S9 : Differences in Matched DFEs ####

# for each strain, get its matched, closest in GR strain in a different env

getClosestStrain_DeltaDFE_adjGR = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice )
  focalStrainGR = focalStrainInfo$adjGR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain !=StrainChoice ) # or to itself in diff env
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }
  
  allPossibleComps$deltaGR = abs(allPossibleComps$adjGR - focalStrainGR )
  
  compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
  
  
  
  ### ~~~ for doing KS test ~~~ ##
  dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
  dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
  
 ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
  
 joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
 
 ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
 
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
  
  
  returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                        deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                        deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                        deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                        KS_p_unjoined = ks_unjoined$statistic,
                        KS_p_joined = ks_joined$statistic)
  
  
  
  
  
  return(returnDF)
}



getClosestStrain_DeltaDFE_rawGR = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice)
  focalStrainGR = focalStrainInfo$Mean_GR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain !=StrainChoice ) # or to itself in diff env
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }
  
  allPossibleComps$deltaGR = abs(allPossibleComps$Mean_GR - focalStrainGR )
  
  compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
 
   ### ~~~ for doing KS test ~~~ ##
  dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
  dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
  
  ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
  
  joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
  
  ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
  
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
  
  
  returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                        deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                        deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                        deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                        KS_p_unjoined = ks_unjoined$statistic,
                        KS_p_joined = ks_joined$statistic)
  
  return(returnDF)
}


getClosestStrain_DeltaDFE_sameStrainDiffEnv = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice)
  focalStrainGR = focalStrainInfo$Mean_GR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain ==StrainChoice ) # ocomp to itself
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }else{
    allPossibleComps$deltaGR = abs(allPossibleComps$Mean_GR - focalStrainGR )
    
    compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
    ### ~~~ for doing KS test ~~~ ##
    dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
    dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
    
    ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
    
    joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
    
    ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
    
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
    
    
    returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                          deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                          deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                          deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                          KS_p_unjoined = ks_unjoined$statistic,
                          KS_p_joined = ks_joined$statistic)
    
    return(returnDF)
  }
 
}



### ~~~ ALL ENVS ##
dfeDatComps_0 = data.frame(Strain = oldDFEdat$Strain, env = oldDFEdat$env)
dfeDatComps_0 = subset(dfeDatComps_0, env != 'YPD')

dfeDatComps_adjGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_adjGR(Strain,env)))%>%
  unnest_wider(output)
  
dfeDatComps_adjGR$type = 'adjGR'

dfeDatComps_rawGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_rawGR(Strain,env)))%>%
  unnest_wider(output)

dfeDatComps_rawGR$type = 'rawGR'


dfeDatComps_sameStrain = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_sameStrainDiffEnv(Strain,env)))%>%
  unnest_wider(output)

dfeDatComps_sameStrain$type = 'sameStrainDiffEnv'

joined_dfeComp = rbind(dfeDatComps_adjGR,dfeDatComps_rawGR,dfeDatComps_sameStrain )


# Now, look at stats of DFE sim 
meandf = joined_dfeComp %>% 
  group_by(type) %>% 
  summarize(meanVal_u = mean(deltaMean), 
            meanSD = mean(deltaVar), 
            meanSkew = mean(deltaSkew),
            meanKS = mean(KS_p_unjoined))



pMean = ggplot(joined_dfeComp, aes(x = deltaMean, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanVal_u, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,80), breaks = c(0,40,80))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


pVar = ggplot(joined_dfeComp, aes(x = deltaVar, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSD, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,80), breaks = c(0,40,80))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)



pSkew = ggplot(joined_dfeComp, aes(x = deltaSkew, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSkew, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,80), breaks = c(0,40,80))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


pKS = ggplot(joined_dfeComp, aes(x = KS_p_unjoined, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanKS, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,80), breaks = c(0,40,80))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


comp_DFEs_meanGR_vs_adjGR_allenv = ggarrange(plotlist = list(pMean,pVar,pSkew,pKS), nrow = 1, ncol = 4)

#ggsave(paste0(figSave,'comp_DFEs_meanGR_vs_adjGR_allEnvs.pdf' ),comp_DFEs_meanGR_vs_adjGR_allenv, width = 17, height = 5, unit = 'cm')





### ~~~ Just most different envs ##
dfeDatComps_0 = data.frame(Strain = oldDFEdat$Strain, env = oldDFEdat$env)
dfeDatComps_0 = subset(dfeDatComps_0, env %in% c('30SC5', '37SC7'))

dfeDatComps_adjGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_adjGR(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_adjGR$type = 'adjGR'

dfeDatComps_rawGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_rawGR(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_rawGR$type = 'rawGR'


dfeDatComps_sameStrain = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_sameStrainDiffEnv(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_sameStrain$type = 'sameStrainDiffEnv'

joined_dfeComp = rbind(dfeDatComps_adjGR,dfeDatComps_rawGR,dfeDatComps_sameStrain )



# Now, look at stats of DFE sim 
meandf = joined_dfeComp %>% 
  group_by(type) %>% 
  summarize(meanVal_u = mean(deltaMean, na.rm = T), 
            meanSD = mean(deltaVar, na.rm = T), 
            meanSkew = mean(deltaSkew, na.rm = T),
            meanKS = mean(KS_p_unjoined, na.rm = T))



pMean = ggplot(joined_dfeComp, aes(x = deltaMean, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanVal_u, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(0,15,30))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


pVar = ggplot(joined_dfeComp, aes(x = deltaVar, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSD, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(0,15,30))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)



pSkew = ggplot(joined_dfeComp, aes(x = deltaSkew, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSkew, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(0,15,30))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


pKS = ggplot(joined_dfeComp, aes(x = KS_p_unjoined, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanKS, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(0,15,30))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.title.y =  element_blank(), 
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_blank(), 
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

comp_DFEs_meanGR_vs_adjGR = ggarrange(plotlist = list(pMean,pVar,pSkew,pKS), nrow = 1, ncol = 4)

#ggsave(paste0(figSave,'comp_DFEs_meanGR_vs_adjGR_envsMostDiff.pdf' ),comp_DFEs_meanGR_vs_adjGR, width = 17, height = 5, unit = 'cm')







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S10 : Effect of Missing Measures on DFE moments ####
# ~~ Here, went through each environment, clustered mutations for whether have a measurement in each strain
# and then manually chose a set of 40 mutations measured in as many strains as possible

## ~~~~ 30SC3 
#paste(shQuote(levels(orderedDF$Mut_ID)[28:67]), collapse=", ") # get 40 mutations
# c('3', '121', '59', '123', '115', '105', '67', '68', '1', '85', '103', '84', '55', '127', '124', '114', '110', '108', '104', '101', '92', '81', '80', '72', '66', '58', '44', '35', '17', '15', '14', '10', '2', '7', '18', '46', '11', '62', '29', '88')

SubStrainEnv_30SC3 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC3' & Mut_ID %in% c('3', '121', '59', '123', '115', '105', '67', '68', '1', '85', '103', '84', '55', '127', '124', '114', '110', '108', '104', '101', '92', '81', '80', '72', '66', '58', '44', '35', '17', '15', '14', '10', '2', '7', '18', '46', '11', '62', '29', '88'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC3 = SubStrainEnv_30SC3 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC3 = subset(SubStrainEnv_30SC3, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC3$Strain)) # get 30 strains
#write.csv(SubStrainEnv_30SC3, paste0(figSave, 'SubStrainEnv_30SC3.csv'))


## ~~~~ 30SC5 
#paste(shQuote(levels(orderedDF$Mut_ID)[25:64]), collapse=", ") # get 40 mutations
# c('18', '127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '92', '88', '82', '81', '80', '72', '70', '68', '66', '62', '58', '55', '46', '44', '35', '29', '27', '17', '15', '14', '11', '10', '7', '2', '3', '97', '21', '22', '69')

SubStrainEnv_30SC5 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC5' & Mut_ID %in%  c('18', '127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '92', '88', '82', '81', '80', '72', '70', '68', '66', '62', '58', '55', '46', '44', '35', '29', '27', '17', '15', '14', '11', '10', '7', '2', '3', '97', '21', '22', '69'))

# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC5 = SubStrainEnv_30SC5 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC5 = subset(SubStrainEnv_30SC5, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC5$Strain)) # get 36 strains
#write.csv(SubStrainEnv_30SC5, paste0(figSave, 'SubStrainEnv_30SC5.csv'))


## ~~~~ 30SC7 
#paste(shQuote(levels(orderedDF$Mut_ID)[43:82]), collapse=", ") # get 40 mutations
# c('127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '103', '101', '92', '88', '85', '82', '81', '80', '72', '70', '66', '62', '58', '55', '46', '44', '35', '29', '27', '18', '17', '14', '11', '10', '2', '7', '69', '3', '15', '22')

SubStrainEnv_30SC7 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC7' & Mut_ID %in%  c('127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '103', '101', '92', '88', '85', '82', '81', '80', '72', '70', '66', '62', '58', '55', '46', '44', '35', '29', '27', '18', '17', '14', '11', '10', '2', '7', '69', '3', '15', '22'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC7 = SubStrainEnv_30SC7 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC7 = subset(SubStrainEnv_30SC7, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC7$Strain)) # get 40 strains
#write.csv(SubStrainEnv_30SC7, paste0(figSave, 'SubStrainEnv_30SC7.csv'))


## ~~~~ 37SC3 
#paste(shQuote(levels(orderedDF$Mut_ID)[23:62]), collapse=", ") # get 40 mutations
# c('27', '85', '105', '11', '26', '82', '121', '115', '68', '88', '92', '122', '17', '55', '84', '110', '15', '80', '1', '101', '83', '46', '18', '127', '124', '114', '108', '81', '72', '62', '44', '35', '29', '2', '10', '7', '104', '123', '50', '3')

SubStrainEnv_37SC3 = subset(df_mutsCanFitLine_wAdjGR, env == '37SC3' & Mut_ID %in%  c('27', '85', '105', '11', '26', '82', '121', '115', '68', '88', '92', '122', '17', '55', '84', '110', '15', '80', '1', '101', '83', '46', '18', '127', '124', '114', '108', '81', '72', '62', '44', '35', '29', '2', '10', '7', '104', '123', '50', '3'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC3 = SubStrainEnv_37SC3 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_37SC3 = subset(SubStrainEnv_37SC3, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC3$Strain)) # get 23 strains
#write.csv(SubStrainEnv_37SC3, paste0(figSave, 'SubStrainEnv_37SC3.csv'))


## ~~~~ 37SC5 
#paste(shQuote(levels(orderedDF$Mut_ID)[8:47]), collapse=", ") # get 40 mutations
# c('92', '82', '105', '84', '26', '55', '67', '123', '122', '83', '50', '58', '127', '27', '115', '18', '121', '1', '104', '101', '7', '85', '72', '81', '62', '80', '46', '110', '68', '15', '124', '114', '108', '44', '35', '29', '17', '2', '10', '60')

SubStrainEnv_37SC5 = subset(df_mutsCanFitLine_wAdjGR, env == '37SC5' & Mut_ID %in%  c('92', '82', '105', '84', '26', '55', '67', '123', '122', '83', '50', '58', '127', '27', '115', '18', '121', '1', '104', '101', '7', '85', '72', '81', '62', '80', '46', '110', '68', '15', '124', '114', '108', '44', '35', '29', '17', '2', '10', '60'))

# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC5 = SubStrainEnv_37SC5 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_37SC5 = subset(SubStrainEnv_37SC5, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC5$Strain)) # get 16 strains
#write.csv(SubStrainEnv_37SC5, paste0(figSave, 'SubStrainEnv_37SC5.csv'))




### ~~~~~ 37SC7 
# c('LK4-B12', 'LK2-E08', 'LK2-D05', 'LK4-B01', 'LK5-F08', 'LK5-H12', 'LK3-C04', 'LK2-B07', 'LK3-B08', 'LK3-C12', 'LK5-D09', 'LK2-D07', 'LK6-A05', 'LK1-C09', 'LK4-H11', 'LK2-A10', 'LK5-G01', 'LK1-B05', 'LK5-G03', 'LK2-A06', 'LK3-H11')
# c('82', '88', '121', '123', '84', '58', '92', '50', '27', '29', '3', '55', '1', '11', '83', '105', '62', '72', '18', '104', '101', '114', '35', '124', '108', '81', '44', '17', '15', '10', '2', '7', '26', '127', '70', '80', '85', '46', '122', '103')

#paste(shQuote(levels(orderedDF$name)[1:21]), collapse=", ")
#paste(shQuote(levels(orderedDF$Mut_ID)[1:50]), collapse=", ")

SubStrainEnv_37SC7 = subset(df_mutsCanFitLine_wAdjGR,env == '37SC7' &  Mut_ID %in% c('82', '88', '121', '123', '84', '58', '92', '50', '27', '29', '3', '55', '1', '11', '83', '105', '62', '72', '18', '104', '101', '114', '35', '124', '108', '81', '44', '17', '15', '10', '2', '7', '26', '127', '70', '80', '85', '46', '122', '103'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC7 = SubStrainEnv_37SC7 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))
SubStrainEnv_37SC7 = subset(SubStrainEnv_37SC7, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC7$Strain)) # get 16 strains
#write.csv(SubStrainEnv_37SC7, paste0(figSave, 'SubStrainEnv_37SC7.csv'))




### ~~~~~~ PUT ALL TOGETHER

allSubDat = rbind(SubStrainEnv_30SC3,SubStrainEnv_30SC5,SubStrainEnv_30SC7,SubStrainEnv_37SC3,SubStrainEnv_37SC5,SubStrainEnv_37SC7)
noMissingMuts_s_ests = allSubDat

dfeQualCheck_DF = allSubDat %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err, F0) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T))

dfeQualCheck_DF$predDFEmean = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEmean
dfeQualCheck_DF$predDFEvar = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEvar
dfeQualCheck_DF$predDFEskew = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEskew


dfeMean_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ),shape = 16, color = 'white')+
  geom_point(aes(y = DFEmean ),shape = 16, alpha = 0.5)+
  geom_smooth(aes(y = DFEmean ), method = 'lm',formula = 'y~x', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




dfeVar_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR ,color = env))+
  #  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ),shape = 16, color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), shape = 16,alpha = 0.5)+
  geom_smooth(aes(y = DFEvar ), method = 'lm',formula = 'y~poly(x,2)', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )






dfeQualCheck_DF$sqGR = dfeQualCheck_DF$adjGR^2
lmVar = summary(lm(DFEvar~adjGR+sqGR, data =dfeQualCheck_DF ))
DFEvar_p0 = lmVar$coefficients[1,1] * (lmVar$coefficients[1,4]<0.05)
DFEvar_p1 =  lmVar$coefficients[2,1] * (lmVar$coefficients[2,4]<0.05)
DFEvar_p2 = lmVar$coefficients[3,1] * (lmVar$coefficients[3,4]<0.05)
DFEvar_p0_se =  lmVar$coefficients[1,2] 
DFEvar_p1_se =  lmVar$coefficients[2,2]
DFEvar_p2_se =  lmVar$coefficients[3,2]



dfeQualCheck_DF$newPredVar = DFEvar_p0 + DFEvar_p1*dfeQualCheck_DF$adjGR + DFEvar_p2*dfeQualCheck_DF$sqGR
dfeQualCheck_DF$cuGR = dfeQualCheck_DF$adjGR^3
dfeQualCheck_DF$cuGR_scaled = dfeQualCheck_DF$cuGR / (dfeQualCheck_DF$newPredVar + dfeQualCheck_DF$sqGR)^(3/2) ##?? should I scal by actual gr or pred?? I think pred bc is already fit
dfeQualCheck_DF$sqGR_scaled = dfeQualCheck_DF$sqGR / (dfeQualCheck_DF$newPredVar + dfeQualCheck_DF$sqGR)^(3/2) 
dfeQualCheck_DF$GR_scaled = dfeQualCheck_DF$adjGR /(dfeQualCheck_DF$newPredVar+ dfeQualCheck_DF$sqGR)^(3/2) 


lmSkew = summary(lm(DFEskew~GR_scaled+sqGR_scaled+cuGR_scaled, data =dfeQualCheck_DF ))
DFE_3_p0 = lmSkew$coefficients[1,1] * (lmSkew$coefficients[1,4]<0.05)
DFE_3_p1 =  lmSkew$coefficients[2,1] * (lmSkew$coefficients[2,4]<0.05)
DFE_3_p2 = lmSkew$coefficients[3,1] * (lmSkew$coefficients[3,4]<0.05)
DFE_3_p3 = lmSkew$coefficients[4,1] * (lmSkew$coefficients[4,4]<0.05)
DFE_3_p0_se =  lmSkew$coefficients[1,2] 
DFE_3_p1_se =  lmSkew$coefficients[2,2]
DFE_3_p2_se =  lmSkew$coefficients[3,2]
DFE_3_p3_se =  lmSkew$coefficients[4,2]

dfeQualCheck_DF$predSkew_new = DFE_3_p0 + DFE_3_p1*dfeQualCheck_DF$GR_scaled + DFE_3_p2*dfeQualCheck_DF$sqGR_scaled + DFE_3_p3*dfeQualCheck_DF$cuGR_scaled


dfeSkew_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ), shape = 16,color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEskew ),shape = 16, alpha = 0.5)+
#  geom_smooth(aes(y = DFEskew ), method = 'lm',formula = 'y~poly(x,3)', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  geom_line(aes(y = predSkew_new), color = 'black', linewidth = 1,lineend='round', linetype = 'dashed')+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )



# ggsave(paste0(figSave,'dfeMean_qualCheck.pdf' ),dfeMean_qualCheck, width = 5 , height =4, unit = 'cm')
# ggsave(paste0(figSave,'dfeVar_qualCheck.pdf' ),dfeVar_qualCheck, width = 5, height = 4, unit = 'cm')
# ggsave(paste0(figSave,'dfeSkew_qualCheck.pdf' ),dfeSkew_qualCheck, width = 5, height = 4, unit = 'cm')
# 



# ~~~ Number of Measurements per Strain cor with GR #####
checkNumMutsMeasured = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR) %>%
  summarize(numMutsMeasured = length(unique(Mut_ID)))


plotNumMutPerStrain = ggplot(checkNumMutsMeasured, aes(x = Mean_GR, y =numMutsMeasured, color = env ))+
  geom_point(shape = 16, alpha = 0.5)+
  geom_smooth(method = 'lm', se = F)+
  xlab("Background GR")+
  ylab('Number of Measured Mutations')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        legend.position = 'none' )
# ggsave(paste0(figSave,'plotNumMutPerStrain.pdf' ),plotNumMutPerStrain, width = 3, height = 3, unit = 'cm')



checkNumMutsMeasured %>%
  nest(data = -(env))%>%
  mutate(lmout = map(data, ~lm(numMutsMeasured~Mean_GR, .)), glance = map(lmout, ~glance(.))) %>%
  unnest(glance)

# ~~~ Hist of number of mutations per strain ####
dfe_dat_test = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T), nPoints = sum(!is.na(avgS)))


numMutsPerStrain = ggplot(dfe_dat_test, aes(x = nPoints))+
  stat_binline(bins = 9, fill = 'transparent', aes(y = after_stat(count)))+
  theme_classic()+
  geom_vline(xintercept = median(dfe_dat_test$nPoints), linewidth = 0.5, color = 'grey')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        legend.position = 'none' )

#ggsave(paste0(figSave,'numMutsPerStrain.pdf' ),numMutsPerStrain, width = 3, height = 3, unit = 'cm')


# ~~~ Reduced Strain Set ####
#remove all strains with < 60 muts 


dfe_dat_test = subset(dfe_dat_test, nPoints >= 60)

dfe_dat_test$predDFEmean = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEmean
dfe_dat_test$predDFEvar = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEvar
dfe_dat_test$predDFEskew = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEskew



p1 =ggplot(dfe_dat_test, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ),shape = 16, color = 'white')+
  geom_point(aes(y = DFEmean ),shape = 16, alpha = 0.5)+
   geom_smooth(aes(y =DFEmean),color = 'black',method = 'lm', se = F, linetype = 'dashed')+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )



p2 = ggplot(dfe_dat_test, aes(x = adjGR ,color = env))+
  #  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ), shape = 16,color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), shape = 16,alpha = 0.5)+
  geom_smooth(aes(y = DFEvar ), method = 'lm',formula = 'y~poly(x,2)', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




dfe_dat_test$sqGR = dfe_dat_test$adjGR^2
lmVar = summary(lm(DFEvar~adjGR+sqGR, data =dfe_dat_test ))
DFEvar_p0 = lmVar$coefficients[1,1] * (lmVar$coefficients[1,4]<0.05)
DFEvar_p1 =  lmVar$coefficients[2,1] * (lmVar$coefficients[2,4]<0.05)
DFEvar_p2 = lmVar$coefficients[3,1] * (lmVar$coefficients[3,4]<0.05)
DFEvar_p0_se =  lmVar$coefficients[1,2] 
DFEvar_p1_se =  lmVar$coefficients[2,2]
DFEvar_p2_se =  lmVar$coefficients[3,2]



dfe_dat_test$newPredVar = DFEvar_p0 + DFEvar_p1*dfe_dat_test$adjGR + DFEvar_p2*dfe_dat_test$sqGR
dfe_dat_test$cuGR = dfe_dat_test$adjGR^3
dfe_dat_test$cuGR_scaled = dfe_dat_test$cuGR / (dfe_dat_test$newPredVar + dfe_dat_test$sqGR)^(3/2) ##?? should I scal by actual gr or pred?? I think pred bc is already fit
dfe_dat_test$sqGR_scaled = dfe_dat_test$sqGR / (dfe_dat_test$newPredVar + dfe_dat_test$sqGR)^(3/2) 
dfe_dat_test$GR_scaled = dfe_dat_test$adjGR /(dfe_dat_test$newPredVar+ dfe_dat_test$sqGR)^(3/2) 


lmSkew = summary(lm(DFEskew~GR_scaled+sqGR_scaled+cuGR_scaled, data =dfe_dat_test ))
DFE_3_p0 = lmSkew$coefficients[1,1] * (lmSkew$coefficients[1,4]<0.05)
DFE_3_p1 =  lmSkew$coefficients[2,1] * (lmSkew$coefficients[2,4]<0.05)
DFE_3_p2 = lmSkew$coefficients[3,1] * (lmSkew$coefficients[3,4]<0.05)
DFE_3_p3 = lmSkew$coefficients[4,1] * (lmSkew$coefficients[4,4]<0.05)
DFE_3_p0_se =  lmSkew$coefficients[1,2] 
DFE_3_p1_se =  lmSkew$coefficients[2,2]
DFE_3_p2_se =  lmSkew$coefficients[3,2]
DFE_3_p3_se =  lmSkew$coefficients[4,2]
dfe_dat_test$predSkew_new = DFE_3_p0 + DFE_3_p1*dfe_dat_test$GR_scaled + DFE_3_p2*dfe_dat_test$sqGR_scaled + DFE_3_p3*dfe_dat_test$cuGR_scaled


p3 =ggplot(dfe_dat_test, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ),shape = 16, color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEskew ),shape = 16, alpha = 0.5)+
#  geom_smooth(aes(y = DFEskew ), method = 'lm',formula = 'y~poly(x,3)', se = F, linetype = 'dashed',lineend='round', color = 'black')+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predSkew_new), color = 'black', linewidth = 1,lineend='round', linetype = 'dashed')+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )

#ggarrange(plotlist = list(p1,p2,p3), nrow =1, ncol = 3)

# ggsave(paste0(figSave,'dfeMean_qualCheck_reducedStrains.pdf' ),p1,  width = 5 , height =4, unit = 'cm')
# ggsave(paste0(figSave,'dfeVar_qualCheck_reducedStrains.pdf' ),p2,  width = 5 , height =4, unit = 'cm')
# ggsave(paste0(figSave,'dfeSkew_qualCheck_reducedStrains.pdf' ),p3,  width = 5 , height =4, unit = 'cm')




# ~~~ Imputed mutation effects ####
allMut_Env_Strain_Combos = expand.grid(Mut_ID = unique(df_s_ests_wGR$Mut_ID),
                                       env = unique(df_s_ests_wGR$env),
                                       Strain = unique(df_s_ests_wGR$Strain))

df_s_ests_wGR_wMissingMuts = full_join(allMut_Env_Strain_Combos,df_s_ests_wGR, by = c('Strain', 'env', 'Mut_ID') )

# add EQ(2) slopes and ints to the df

df_s_ests_wGR_wMissingMuts = full_join(df_s_ests_wGR_wMissingMuts, df_Mutslopes_1slope_6int, by = c('Mut_ID', 'env'))

# add in actual GR dat
df_s_ests_wGR_wMissingMuts = df_s_ests_wGR_wMissingMuts[, !(names(df_s_ests_wGR_wMissingMuts) == 'Mean_GR')]
df_s_ests_wGR_wMissingMuts = full_join(df_s_ests_wGR_wMissingMuts, growthRate_data, by = c('Strain', 'env'))

just37missing = subset(df_s_ests_wGR_wMissingMuts, env %in% c('37SC3', '37SC5', '37SC7') & is.na(avgS) )
dim(just37missing) # 4252 missing
dim(subset(just37missing, is.na(slope))) # 684 of these dont have slope, so are not imputatble
dim(subset(just37missing, !is.na(slope))) # 3568 are  imputatble


## can't impute for the ones I dont have a slope and int for
df_s_ests_wGR_wMissingMuts = subset(df_s_ests_wGR_wMissingMuts, !is.na(slope) & !is.na(intercept) & !is.na(Mean_GR))
# removes 1227 rows of the 6634 missing measures (18%)


df_s_ests_wGR_wMissingMuts$wasMissing = is.na(df_s_ests_wGR_wMissingMuts$avgS) # indicator for later 
df_s_ests_wGR_wMissingMuts_haveS = subset(df_s_ests_wGR_wMissingMuts, !is.na(avgS))
df_s_ests_wGR_wMissingMuts_noS = subset(df_s_ests_wGR_wMissingMuts, is.na(avgS))


## just mean predicted s at all point,, no res var

test = df_s_ests_wGR_wMissingMuts_noS
test$avgS = test$Mean_GR*test$slope + test$intercept

test2 = rbind(df_s_ests_wGR_wMissingMuts_haveS, test)
subdat = subset(test2, env %in% c('37SC3', '37SC5', '37SC7'))
meanMissing = mean(subset(subdat, wasMissing == T)$avgS)
meanNonMissing = mean(subset(subdat, wasMissing == F)$avgS)
missingFitEffect_pooled = ggplot(subset(test2, env %in% c('37SC3', '37SC5', '37SC7')), aes(x = avgS, color = wasMissing))+
  stat_binline(binwidth = 0.03, aes(color = wasMissing,y = after_stat(density)),geom="step", position = position_dodge(width = 0.001),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('darkgreen', 'purple'), drop = FALSE)+
  geom_vline(xintercept = meanMissing, color = 'purple')+
  geom_vline(xintercept = meanNonMissing, color = 'darkgreen')+
  #  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96' ,'#6a9c67'), drop = FALSE)+
  # facet_wrap(~env)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.55, lineend = 'round')+
  theme_classic()+
  xlab('Fitness Effect')+
  ylab('Density')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.position = 'none')

# ggsave(paste0(figSave, 'missingFitEffect_pooled.pdf'), missingFitEffect_pooled , unit = 'cm', width = 3.7, height = 4)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# Figure S11 : QTL analysis ####

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


# ~~ Find Bloom QTLs ####


allGeneDiffs_wide$ncharrow1 = nchar(allGeneDiffs_wide$segregant.27915_chr01_27915_T_C)
allGeneDiffs_wide$X27915_chr01_27915_T_C = substr(allGeneDiffs_wide$segregant.27915_chr01_27915_T_C, allGeneDiffs_wide$ncharrow1,allGeneDiffs_wide$ncharrow1 )
allGeneDiffs_wide = allGeneDiffs_wide[, !(names(allGeneDiffs_wide)== 'ncharrow1')]


# ~~ Clean bloom data
names(bloom_temp_pH_qtls) = c('Trait','Frac_Phen_Var_exp', 'Chromosome', 'Peak_position', 'LOD', 'Left_CI', 'Right_CI', 'Genes_in_CI')
bloom_temp_pH_qtls$csome_alt = bloom_temp_pH_qtls$Chromosome
bloom_temp_pH_qtls$csome_alt[bloom_temp_pH_qtls$Chromosome < 10] = paste0('0', bloom_temp_pH_qtls$Chromosome[bloom_temp_pH_qtls$Chromosome < 10])
bloom_temp_pH_qtls$qtl = paste0('X', bloom_temp_pH_qtls$Right_CI, '_chr', bloom_temp_pH_qtls$csome_alt, '_', bloom_temp_pH_qtls$Right_CI)
bloom_temp_pH_qtls = subset(bloom_temp_pH_qtls, !(Trait %in% c('YNB',"YPD:15C", "YPD:4C" )))# only keep pH or extreme high temp qtls

# ~~ pick top 2 QTls per env
bloom_temp_pH_qtls$qtl_alt = paste0('chr',bloom_temp_pH_qtls$Chromosome, '_', bloom_temp_pH_qtls$Peak_position)

bloom_temp_pH_qtls_sub = bloom_temp_pH_qtls %>%
  group_by(Trait) %>%
  mutate(Rank = rank(as.numeric(-1*Frac_Phen_Var_exp), ties.method = 'first')) %>%# -1 so does decreasing rank
  filter(Rank %in% c(1:2)) 



# ~~ make gene diffs long
allGeneDiffs = gather(allGeneDiffs_wide, key = 'qtl_long', value = 'val', "X28652_chr01_28652_G_T":"X27915_chr01_27915_T_C")
allGeneDiffs$nchar = nchar(allGeneDiffs$qtl_long)
allGeneDiffs$qtl = substr(allGeneDiffs$qtl_long, 1, (allGeneDiffs$nchar - 4))
allGeneDiffs$position = sub(".*_(\\d+)$", "\\1", allGeneDiffs$qtl)
## ~~ need to fix ones that don't work right (have extra _ later)
allGeneDiffs$position = as.numeric(allGeneDiffs$position) # convert ones that dint work to NA
allGeneDiffs$position[is.na(allGeneDiffs$position)] =   sub(".*_(\\d+)_.*$", "\\1", allGeneDiffs$qtl[is.na(allGeneDiffs$position)])


allGeneDiffs$chromosome =sub(".*chr(\\d+)_.*", "\\1", allGeneDiffs$qtl)
allGeneDiffs$chromosome = as.numeric(allGeneDiffs$chromosome)
allGeneDiffs$rowID = 1:nrow(allGeneDiffs)
allGeneDiffs$segregant.27915_chr01_27915_T_C = substr(allGeneDiffs$segregant.27915_chr01_27915_T_C, 1,7)
names(allGeneDiffs)[1] = 'segregant'


## Find locus closest to QTL peak
qtls = c('KRE33_chr14_376315','chr14_470303','chr12_646707','chr15_154799')
# also in this data    # close, but tech diff # close but tech diff # VERY close, but tech diff

allLociChoice = data.frame()
for(i in 1:nrow(bloom_temp_pH_qtls_sub))
{
  row = bloom_temp_pH_qtls_sub[i, ]
  subGenes = subset(allGeneDiffs,chromosome == row$Chromosome )
  subGenes$distFromPeak = abs(as.numeric(subGenes$position) - row$Peak_position)
  locusChoice = subGenes[which.min(subGenes$distFromPeak),]
  
  dfToAdd = locusChoice[, c('qtl', 'position', 'chromosome',  'distFromPeak')]
  dfToAdd$FracVar_bloom = row$Frac_Phen_Var_exp
  dfToAdd$Trait_bloom = row$Trait
  
  
  allLociChoice = rbind(allLociChoice,dfToAdd )
  
}


## CHECK: are any of these significantly correlated with the trait value
# on the same c'some (including sites was using before)
allGeneDiffs$qtl_alt = paste0('chr',allGeneDiffs$chromosome, '_', allGeneDiffs$position)
extraQTLS =  c('chr14_376315','chr14_470303','chr12_646707','chr15_154799',  paste0('chr',allLociChoice$chromosome, '_', allLociChoice$position))

suballGeneDiffs = subset(allGeneDiffs,qtl_alt%in% extraQTLS )

## go thru all pairs
qtlPairs = t(combn(extraQTLS,2)) # all UNIQUE pairs (expand grid gives all pairs when order matters)

allRsq = c()
allP = c()
for(pair in 1:nrow(qtlPairs))
{
  
  sub1 = subset(allGeneDiffs, qtl_alt == qtlPairs[pair,1])
  sub2 = subset(allGeneDiffs, qtl_alt == qtlPairs[pair,2])
  joined = inner_join(sub1, sub2, by = 'segregant')
  
  joined_sum = joined %>%
    group_by(val.x, val.y) %>%
    summarize(num = length(val.x), .groups = 'keep')
  
  # ggplot(joined_sum, aes(x = val.x, y = val.y, fill = num))+
  #   geom_tile()+
  #   geom_text(aes(label=num), color = 'white', size = 10)+
  #   scale_fill_viridis()
  
  out = summary(lm(val.y~val.x, data = joined)) 
  # Rsq = 0.88
  allRsq = c(allRsq, out$adj.r.squared)
  allP = c(allP, out$coefficients[2,4])
}
qtlPairs = as.data.frame(qtlPairs)
qtlPairs$rsq_valSim = allRsq
qtlPairs$p = allP
qtlPairs$adj_p = p.adjust(qtlPairs$p )

subset(qtlPairs, adj_p<0.05 )

## the KRE33 is exactly the same
## the other chromosome 14 one is sig corr with KRE33 (from orig data, so might want to exclude)
## the chromosome 12 ones are sig corr with one another
## all chromosome 15 ones are sig corr,, esp chr15_154799 with chr15_163253 

## So, going to only keep bloom ones , and will only keep loci on diff chromosomes
# chr15_559460: is picked rep from all 15's since is associated with temp and would like to have 1 extra temp associated QTL
# chr14_376315 (KRE33) is the rep from c'some 1 to use
# chr13_769126 : chrom 13 not in old data, so use
# chr12_571414: use the one 12 from the bloom data 

allChosenQTLS = c('chr15_559460','chr14_376315','chr13_769126','chr12_571414')
subset(bloom_temp_pH_qtls_sub, qtl_alt %in% allChosenQTLS)


## So, pull more choices (til can get to 10, on diff c'somes)
bloom_temp_pH_qtls_sub2 = bloom_temp_pH_qtls %>%
  group_by(Trait) %>%
  mutate(Rank = rank(as.numeric(-1*Frac_Phen_Var_exp), ties.method = 'first')) %>%# -1 so does decreasing rank
  filter(Rank %in% c(1:5)) %>%
  filter(!(Chromosome %in% c(12,13,14,15))) %>%
  filter(Peak_position != 385471) #   have decided 4 from before and this has 6 unique chromosomes represented
# pick the one on 4 from from ph8 since have 3 in ph3 and only 2 in 8


allChosenQTLS = c(allChosenQTLS, unique(bloom_temp_pH_qtls_sub2$qtl_alt))


allChosenQTLs_info = bloom_temp_pH_qtls %>%
  filter(qtl_alt %in% allChosenQTLS)%>%
  dplyr::select(Trait, Frac_Phen_Var_exp ,Chromosome ,Peak_position)



## Now, re-find the locus in my gene diffs data associated with each chosen locus

allLociChoice = data.frame()
for(i in 1:nrow(allChosenQTLs_info))
{
  row = allChosenQTLs_info[i, ]
  subGenes = subset(allGeneDiffs,chromosome == row$Chromosome )
  subGenes$distFromPeak = abs(as.numeric(subGenes$position) - row$Peak_position)
  locusChoice = subGenes[which.min(subGenes$distFromPeak),]
  
  dfToAdd = locusChoice[, c('qtl_alt', 'position', 'chromosome',  'distFromPeak')]
  dfToAdd$FracVar_bloom = row$Frac_Phen_Var_exp
  dfToAdd$Trait_bloom = row$Trait
  
  
  allLociChoice = rbind(allLociChoice,dfToAdd )
  
}


## Now, get sub gene differences
suballGeneDiffs = subset(allGeneDiffs,qtl_alt%in% allLociChoice$qtl_alt )
uniqueQTLS = unique( allLociChoice$qtl_alt )

##  ~~~~~ Prop Var exp each locus separately 
names(suballGeneDiffs)[1] ='Strain'
suballGeneDiffs_wGR = inner_join(suballGeneDiffs, growthRate_data, by = 'Strain')

propVar_GR_allLoci_perEnvQTL = function(data)
{
  data_spread = data %>%
    dplyr::select(Mean_GR, qtl_long, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
    pivot_wider(names_from = qtl_long, 
                values_from = val)
  
  jointPropVar_lm = (lm(Mean_GR ~ ., data = data_spread))
  anova_joint =  anova(jointPropVar_lm) 
  afss <- anova_joint$"Sum Sq"
  percExpOut = (cbind(anova_joint,PctExp=afss/sum(afss)*100)) 
  
  jointPropVar = sum(percExpOut$PctExp[rownames(percExpOut)!='Residuals'] * (percExpOut$`Pr(>F)`<0.05)[rownames(percExpOut)!='Residuals'])
  return(data.frame(jointPropVar, rsq =  summary(jointPropVar_lm)$r.squared*(lmp(jointPropVar_lm)<0.05)))
}

propVarExp_eachLocus_eachEnv = suballGeneDiffs_wGR%>%
  nest(data = -c(env,qtl_alt)) %>%
  mutate(pExpAll = map(data, ~propVar_GR_allLoci_perEnvQTL(.)))%>%
  unnest(pExpAll)
# mutate(out = map(data, ~lm(Mean_GR~val,data = .)),  glance = map(out, ~glance(.))) %>%
#unnest(glance) %>%
# dplyr::select(qtl_alt, env, r.squared , adj.r.squared, p.value)

# 
# ggplot(propVarExp_eachLocus_eachEnv, aes(x= env, y= qtl_alt ,fill = jointPropVar))+
#   geom_tile()+
#   geom_text(aes(label=round(jointPropVar)), color = 'white', size = 4)+
#   scale_fill_viridis()

propVarExp_eachLocus_eachEnv$chromosome = as.numeric(sub(".*chr(\\d+)_.*", "\\1", propVarExp_eachLocus_eachEnv$qtl_alt))
propVarExp_eachLocus_eachEnv = propVarExp_eachLocus_eachEnv[order(propVarExp_eachLocus_eachEnv$chromosome), ]
propVarExp_eachLocus_eachEnv$qtl_alt = factor(propVarExp_eachLocus_eachEnv$qtl_alt, levels = unique(propVarExp_eachLocus_eachEnv$qtl_alt))




# ~~~~~ Prop var all loci together per env 
library(car)
propVar_GR_allLoci_perEnv_withFullDat = function(data, envIn )
{
  data_spread = data %>%
    dplyr::select(Mean_GR, qtl_alt, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
    pivot_wider(names_from = qtl_alt, 
                values_from = val)
  
  
  # order the qtls by top to bottom percent var exp
  sub = subset(propVarExp_eachLocus_eachEnv, env == envIn)
  qtlMaxtoMin = as.character(sub$qtl_alt[order(sub$jointPropVar, decreasing = T) ])
  data_spread = data_spread[, c('Mean_GR', qtlMaxtoMin)]
  
  
  jointPropVar_lm = (lm(Mean_GR ~ ., data = data_spread))
  anova_joint =  anova(jointPropVar_lm) 
  afss <- anova_joint$"Sum Sq"
  percExpOut = (cbind(anova_joint,PctExp=afss/sum(afss)*100)) 
  
  jointPropVar = sum(percExpOut$PctExp[rownames(percExpOut)!='Residuals']  * (percExpOut$`Pr(>F)`<0.05)[rownames(percExpOut)!='Residuals'])
  return(data.frame(locus = rownames(percExpOut), p_val = percExpOut$`Pr(>F)`, pctExp = percExpOut$PctExp))
  ##NOTE :: the ANOVA and Rsq results are the SAME tot var explained if I do not remove  N.S predictors
  ## so, I prefer the ANOVA method because removes the spurious effects of N.S factors
  ## more robust than just checking if full model is sig 
  
}

allPropVarJoint_extra = suballGeneDiffs_wGR %>%
  nest(data = -env) %>%
  mutate(envIn = env) %>%
  mutate(pExpAll = map2(data, envIn, ~propVar_GR_allLoci_perEnv_withFullDat(.x, .y))) %>%
  unnest(pExpAll)

allPropVarJoint_extra = subset(allPropVarJoint_extra,locus!='Residuals' )
allPropVarJoint_extra$isSig = allPropVarJoint_extra$p_val<0.05
allPropVarJoint_extra$correctedPctExp = allPropVarJoint_extra$pctExp * allPropVarJoint_extra$isSig


sub = subset(allPropVarJoint_extra, p_val<0.05)
totGRvarExp_colByLocus_alt = ggplot(allPropVarJoint_extra, aes(x = env, y = correctedPctExp, fill = env,color =locus, group = interaction(env, locus)))+
  geom_col()+
  #geom_point(aes(y = propKRE33))+
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  scale_color_manual(values = c(rep('black'), times = 10),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  xlab("Environment")+
  ylab('Total Joint % Variance Explained')+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text.x =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        axis.text.y = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.position = 'none')
#ggsave(paste0(figSave, 'totGRvarExp2.pdf'), totGRvarExp_colByLocus_alt, width = 18, height = 4, unit = 'cm')





testDat_muts_extraQTLs = inner_join(df_mutsCanFitLine_wAdjGR, suballGeneDiffs[,c('Strain', 'val', 'qtl_alt')], by = 'Strain')

## look at all qtls added on mutEffect ####
getp_val_allQTLs_aboveGR_2 = function(data)
{
  data_spread = as.data.frame(data) %>%
    dplyr::select(avgS, Mean_GR,env, qtl_alt, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
    pivot_wider(names_from = qtl_alt, 
                values_from = val)
  # order so test do eq 2 fit first (avgS + GR) then all QTLS
  
  
  # First, assess that there are multiple values at all QTls 
  sdata = data_spread[,c(uniqueQTLS)]
  numUnique_perCol = sapply(sdata, function(x) length(unique(x)))
  
  
  
  
  
  if(sum(numUnique_perCol == 2)< length(numUnique_perCol))
  { # if have no variation in 1 of the qtls
    return(data.frame(sig_GR = NA, P_explained_GR = NA, P_explained_Env = NA,  P_exp_allQTLS_together = NA))
    
  }else{
    
    
    
    # ~~
    lmout = (lm(avgS ~ ., data = data_spread)) 
    af = anova(lmout)
    afss <- af$"Sum Sq"
    percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) 
    
    ## only add signficant ones 
    return(data.frame(sig_GR = percExpOut$`Pr(>F)`[rownames(percExpOut) == 'Mean_GR'], 
                      P_explained_GR = percExpOut$PctExp[rownames(percExpOut) == 'Mean_GR']*(percExpOut$`Pr(>F)`[rownames(percExpOut) == 'Mean_GR']<0.05),  
                      P_explained_Env = (percExpOut$PctExp[rownames(percExpOut) == 'env'])*(percExpOut$`Pr(>F)`[rownames(percExpOut) == 'env']<0.05),
                      P_exp_allQTLS_together = sum(percExpOut$PctExp[!(rownames(percExpOut) %in% c('Mean_GR', 'env', 'Residuals'))] * (percExpOut$`Pr(>F)`[!(rownames(percExpOut) %in% c('Mean_GR', 'env', 'Residuals'))]<0.05 ))  ))
    
    
  }
  
}

allQTls_additive_sig_perMut_env_aboveGR_2 = testDat_muts_extraQTLs %>% 
  nest(data = -c(Mut_ID)) %>%
  mutate(sig_df = map(data, ~getp_val_allQTLs_aboveGR_2(.))[[1]]) %>% 
  unnest(sig_df)

allQTls_additive_sig_perMut_env_aboveGR_2$p_adj_GR = p.adjust(allQTls_additive_sig_perMut_env_aboveGR_2$sig_GR, method = 'BH')
allQTls_additive_sig_perMut_env_aboveGR_2 = allQTls_additive_sig_perMut_env_aboveGR_2[rev(order(allQTls_additive_sig_perMut_env_aboveGR_2$P_exp_allQTLS_together)), ]
allQTls_additive_sig_perMut_env_aboveGR_2$Mut_ID = factor(allQTls_additive_sig_perMut_env_aboveGR_2$Mut_ID, levels = unique(allQTls_additive_sig_perMut_env_aboveGR_2$Mut_ID))

sub = allQTls_additive_sig_perMut_env_aboveGR_2[,c('Mut_ID', 'P_explained_GR','P_explained_Env', 'P_exp_allQTLS_together')]
sub$Mut_ID = as.character(sub$Mut_ID)
## put in mut names
Mut_full_info_Use$Mut_ID = as.factor(Mut_full_info_Use$Mut_ID)
sub = inner_join(sub, Mut_full_info_Use, by = 'Mut_ID')


## put Env+GR = %var explained by Eq2
sub2 = sub
sub2$P_exp_eq2 = sub$P_explained_Env + sub$P_explained_GR
sub2 = sub2[,c('Mut_ID','Gene.Use', 'P_exp_eq2' , 'P_exp_allQTLS_together')]

gathered2 = gather(sub2,key = 'Type', value = 'PctExp',  P_exp_eq2:P_exp_allQTLS_together)
gathered2 = gathered2 %>%
  group_by(Mut_ID) %>%
  mutate(totPctExp = sum(PctExp))
gathered2 = gathered2[rev(order(gathered2$totPctExp)),]

gathered2$Mut_ID = factor(gathered2$Mut_ID, levels = unique(gathered2$Mut_ID))
gathered2$Gene.Use = factor(gathered2$Gene.Use, levels = unique(gathered2$Gene.Use))

gathered2= subset(gathered2, !is.na(totPctExp))
allQTLinfo = gathered2
qtlPlot2 = ggplot(gathered2, aes(x = Gene.Use, y = PctExp, fill = Type))+
  geom_col(color = 'white', linewidth = 0.2)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c( '#c24674', '#46c294'), drop = F)+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5, color = 'black'), 
        axis.text.y = element_text(size = 5),
        legend.position = 'none')
#ggsave(paste0(figSave, 'qtlPlot3.pdf'), qtlPlot2, width = 18, height = 6.4, unit = 'cm')


### pick top 3 loci with % var from qtls
s = subset(gathered2, Type == 'P_exp_allQTLS_together')
s$relProp = s$PctExp/s$totPctExp
t = head(s[order(s$PctExp, decreasing = T), ])
t
# Mut IDs 95,61,57

quantile(subset(s, !(Mut_ID %in% c(95,61,57)))$PctExp)
quantile(subset(s, !(Mut_ID %in% c(95,61,57)))$relProp, na.rm = T)

## 95
data = subset(testDat_muts_extraQTLs, Mut_ID == 95)
data_spread = as.data.frame(data) %>%
  dplyr::select(avgS, Mean_GR,env, qtl_alt, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
  pivot_wider(names_from = qtl_alt, 
              values_from = val)

lmout = (lm(avgS ~., data = data_spread)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100))
percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ] [which.max((percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ])$PctExp), ]
# highest effect qtl is chr4_1006082



## 61
data = subset(testDat_muts_extraQTLs, Mut_ID == 61)
data_spread = as.data.frame(data) %>%
  dplyr::select(avgS, Mean_GR,env, qtl_alt, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
  pivot_wider(names_from = qtl_alt, 
              values_from = val)

lmout = (lm(avgS ~., data = data_spread)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) 
percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ] [which.max((percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ])$PctExp), ]

# highest effect qtl is chr13_769126 


## 57
data = subset(testDat_muts_extraQTLs, Mut_ID == 57)
data_spread = as.data.frame(data) %>%
  dplyr::select(avgS, Mean_GR,env, qtl_alt, val) %>% # use qtlLong instead of qtl_alt so is ok to next with that value
  pivot_wider(names_from = qtl_alt, 
              values_from = val)

lmout = (lm(avgS ~., data = data_spread)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) # env and GR are both NS
percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ] [which.max((percExpOut[!(as.character(rownames(percExpOut)) %in% c('env', "Mean_GR", 'Residuals')), ])$PctExp), ]

# highest effect qtl is chr15_559460


pList = list()
mutsChoose = c( 95,61,57)
corespQTLS = c('chr4_1006082', 'chr13_769126', 'chr15_559460')

for(i in 1:length(mutsChoose))
{
  mChoose = mutsChoose[i]# 80
  #eChoose = sub_sigKre33[i,]$env #'37SC3'
  name = paste0(mChoose, ': ', subset(Mut_full_info_Use, Mut_ID == mChoose )$Gene.Use)
  qtlChoice = corespQTLS[i]
  
  sub = subset(testDat_muts_extraQTLs, Mut_ID == mChoose & qtl_alt ==qtlChoice)
  
  
  # ## get predS
  # best = '1s6i'
  # fitted = resultLRT_Test[resultLRT_Test$Mut_ID == mChoose, paste0('fit_', best)][[1]][[1]]
  # data = subset(df_mutsCanFitLine, Mut_ID == mChoose)
  # data$predS = fitted$fitted.values # note this only works because data here is subsetted from eaxacr same DF in exact same way as when do the LM
  # data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  # 
  fit = lm(avgS ~ Mean_GR + env, data = sub)
  dataFit = data.frame(predS = fit$fitted.values, Mean_GR = sub$Mean_GR, env = sub$env) # can use the meanGR from sub bc in same order fitted
  
  pList[[i]] = ggplot(sub, aes(x = Mean_GR, y = avgS))+
    geom_point(size = 1,  aes(color = env, shape =   as.factor(val), alpha = as.factor(val)))+
    # scale_color_manual(values = c('#f6d53c', '#3c5df6'))+
    geom_line(data = dataFit, aes(x = Mean_GR, y = predS, color = env))+
    scale_color_manual(values = c( '#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'),drop = F)+
    scale_alpha_discrete(range = c(0.34, 0.9))+
    ggtitle(paste0( name, ', ', qtlChoice))+
    geom_hline(yintercept = 0, color = 'grey', linewidth = 1.2)+
    #  ylim(-0.12,0.12)+
    xlim(0,0.4)+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),legend.position = 'none', axis.title = element_blank(),
          axis.text = element_text(color = 'black', size = 8))
}

exampleQTLs = ggarrange(plotlist = pList[c(2,1,3)], nrow = 1, ncol = 3)
#ggsave(paste0(figSave, 'exampleQTLs3.pdf'), exampleQTLs, width= 17, height = 4, unit = 'cm')





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# Figure S12 : 4 site landscape ####



# ~~ Set GR shift and Scale ####
# ~~~ NOTE: these values define how I convert 'fitness' values from the model to realistic fitness values
grScale = 0.2 # usually 0.2 to produce Rsq of slope ~ int = 0.98
grShift = 0.245 # equivalent to the pivot GR, same as 30sc5 in concatenated data, found by running code below then adjusting accordingly
chosenEnv = '30SC5' # choose one environment to do concatenated data regressions with 




# ~~ Define core functions ####
# These functions are exactly as they are implemented in the scripts from 
# Reddy and Desai eLife 2021 https://elifesciences.org/articles/64740
# See the sections on the 'connectedness model' in this paper for more details
# only places they are adjusted is for correct R syntax

# Convert X to n function
convert_X_to_n <- function(X) {
  binvec <- (1 - X) / 2
  binvec <- as.integer(binvec)
  n <- sum(2^(seq_along(X) ) * binvec)
  return(n)
}

# Generate module weights function
generate_module_weights <- function(n_m) {
  return(rnorm(2^n_m)) # sample all epistatic interaction terms
}

# Get module output function
# get all fitness contributions from all interaction levels of these loci
## note: this just implements the Reddy-Desai eLife eq (1)
get_module_output <- function(weights_n, Xn, Ks) {
  Xn = as.matrix(Xn)
  n <- ((1 - (as.matrix(Xn))) / 2) # turn -1,1 to 0,1 vec,, 0 is for 1 and 1 is for -1
  if(ncol(Xn) == 1) # 1 indiv
  {
    coeffs <- (-1)^(Ks %*% as.vector(n)) 
  }else{
    coeffs <- (-1)^(Ks %*% t(n)) # using matrix multiplication and power --- have to have else statement for single individual because the matrix is made in wrong direction with 1 indiv
    
  }
  return(colSums(coeffs * weights_n))
}

# Generate modules function 
generate_modules <- function(mus, M, maxn_m) {
  p <- length(mus)
  weights <- list()
  loci <- list()
  for (m in 1:M) {
    n_m <- 0
    selected_m <- NULL
    while (n_m > maxn_m | n_m < 1) {
      locs <- 1:p
      selected_m <- locs[runif(p) < mus]
      n_m <- length(selected_m)
    }
    weights[[m]] <- generate_module_weights(n_m) / sqrt(M)
    loci[[m]] <- selected_m
  }
  return(list(weights, loci))
}

# Get Ks all function
get_Ks_all <- function(maxn_m) {
  Ks_all <- list()
  for (n_m in 1:maxn_m) {
    Ks <- matrix(0, nrow = 2^n_m, ncol = n_m)
    for (i in 1:(2^n_m )) {
      Ks[i, ] <- as.integer(intToBits(i))[1:n_m] # Convert binary to integer array
    }
    Ks_all[[n_m ]] <- Ks
  }
  return(Ks_all)
}

# Get fitness function
get_fitness <- function(weights, loci, X, Ks_all) {
  fitness <- numeric(nrow(X))
  for (m in 1:M) {
    n_m <- length(loci[[m]])
    Ks <- Ks_all[[n_m ]]
    fitness <- fitness + get_module_output(weights[[m]], X[, loci[[m]]], Ks)
  }
  return(fitness)
}

# Compute Vloci function
compute_Vloci <- function(weights, loci, Ks_all, p) {
  Vloci <- numeric(p)
  V <- 0
  for (m in 1:length(weights)) {
    n_m <- length(loci[[m]])
    Ks <- Ks_all[[n_m + 1]]
    V <- V + sum(weights[[m]][-1]^2)  # Exclude the first element and square the rest
    for (i in 1:n_m) {
      filt <- which(Ks[, i] == 1)
      Vloci[loci[[m]][i]] <- Vloci[loci[[m]][i]] + sum(weights[[m]][filt + 1]^2)
    }
  }
  Vloci <- Vloci / V
  return(Vloci)
}



# Evolve function
evolve_ngens <- function(indiv, nGens, nIts, U, N, weights, loci, Ks_all) {
  evolved <- list()
  gen = num
  Ys <- (get_fitness(weights, loci, indiv, Ks_all) * grScale) +  grShift
  
  allEvo_oneIt = list()
  gensOfFix_oneIt = c()
  
  df_allIt = data.frame()
  
  for(it in 1:nIts)
  {
    print(it)
    times <- 1:nGens
    fixloci <- numeric(nGens)
    Ytraj <- numeric(nGens)+Ys
    vselected <- numeric(nGens)
    diffFromFounder = numeric(nGens)
    numNegOnes = numeric(nGens)+sum(indiv == -1)
    
    ## save all individuals from last iteration
    allIndiv_1it = t(matrix(indiv))
    
    listOfFixedMuts = c()
    countNumRepeatSites = 0
    nFixed = 0
    
    founder <- indiv
    curr <- founder
    Ycurr <- Ys
    
    i = 0
    # Sample waiting time until next mutation 
    waitingTime = ceiling(rexp(1, rate = N*U))
    
    while (i < nGens) 
    {
      i = i + waitingTime
      temp <- curr
      
      # sample locus to mutate 
      locus <- sample(1:p, 1)
      
      # change sign at this locus 
      temp[locus] <- temp[locus] * -1
      
      # calculate fitness of this new genotype 
      Ynew <- get_fitness(weights, loci, matrix(temp, nrow = 1), Ks_all)[1] * grScale + grShift
      s <- Ynew - Ycurr
      if (runif(1) < 2 * s) { # fixation probability =  2s 
        Ycurr <- Ynew
        curr <- temp
        fixloci[i:nGens] <-  locus
        Ytraj[i:nGens] <-  Ycurr
        vselected[i:nGens] <-  Vloci[locus]
        diffFromFounder[i:nGens] = sum(curr != founder) # only need to update when get new fix bc writes same value for all future gens
        numNegOnes[i:nGens] = sum(curr == -1)
        
        
        allIndiv_1it = rbind(allIndiv_1it,curr )
        
        
        countNumRepeatSites = countNumRepeatSites + locus %in% listOfFixedMuts
        listOfFixedMuts = c(listOfFixedMuts, locus)
        nFixed = nFixed+1
        
        if(it == nIts) # save all data for 1 environment 
        {
          allEvo_oneIt[[i]] = curr
          gensOfFix_oneIt = c(gensOfFix_oneIt, i)
        }
        
      }else{ # if it does not fix 
        fixloci[i:nGens] <-  NA
        Ytraj[i:nGens] <-  Ycurr
        vselected[i:nGens] <-  NA
        
        
      }
      
      
      # Sample time until next mutation
      waitingTime =  ceiling(rexp(1, rate = N*U))
      
    }
    
    evolved[[it]] = curr # save final indiviual 
    
    dfIt = data.frame(times, fixloci[1:nGens], Ytraj[1:nGens], vselected[1:nGens], diffFromFounder[1:nGens], countNumRepeatSites, numNegOnes[1:nGens])
    names(dfIt) = c('time', 'fixloci', 'fitness', 'vselected', 'diffFromFounder', 'totalNumRepeatSites', 'numNeg1')
    dfIt$it = it
    
    df_allIt = rbind(df_allIt, dfIt)
  }
  
  
  return(list(df_allIt,evolved,allIndiv_1it, allEvo_oneIt,gensOfFix_oneIt))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# Preliminary analysis of our data ####
# To do this, calculate concatenated data slopes and intercepts 
# Implemented as in the Reddy and Desai description of 
# "Analysis of Data from Johnson et al 2019" section 
# in the materials and methods

# ~~ Step 1: Concat df_s_ests data ####
df_s_ests_wGR_insertions = df_s_ests_wGR
df_s_ests_wGR_insertions$type = 'insertion'

df_s_ests_wGR_reversions = df_s_ests_wGR
df_s_ests_wGR_reversions$Mean_GR = df_s_ests_wGR$Mean_GR + df_s_ests_wGR$avgS
df_s_ests_wGR_reversions$avgS = -1*df_s_ests_wGR$avgS
df_s_ests_wGR_reversions$Std_err = sqrt(df_s_ests_wGR_reversions$Std_err^2 + df_s_ests_wGR_reversions$seS^2) # error prop in the addition for new error in bg GR (can assume same std err for s)
df_s_ests_wGR_reversions$type = 'reversion'


# update mutation calls for the reversions
df_s_ests_wGR_reversions$mutCall = NA
df_s_ests_wGR_reversions$mutCall = NA
df_s_ests_wGR_reversions$minCI = df_s_ests_wGR_reversions$avgS - 2.58*df_s_ests_wGR_reversions$seS ## assuming se of reversion is same as the insertion
df_s_ests_wGR_reversions$maxCI = df_s_ests_wGR_reversions$avgS + 2.58*df_s_ests_wGR_reversions$seS

df_s_ests_wGR_reversions$mutCall[df_s_ests_wGR_reversions$minCI > 0] = 'Beneficial'
df_s_ests_wGR_reversions$mutCall[df_s_ests_wGR_reversions$maxCI < 0] = 'Deleterious'
df_s_ests_wGR_reversions$mutCall[df_s_ests_wGR_reversions$maxCI > 0 & df_s_ests_wGR_reversions$minCI < 0] = 'Neutral'


# concat
allMutEffects = rbind(df_s_ests_wGR_insertions, df_s_ests_wGR_reversions)
remove(df_s_ests_wGR_reversions,df_s_ests_wGR_insertions )

## SUBSET TO FOCAL ENVIRONMENT 
allMutEffects_allenvs = allMutEffects
allMutEffects = subset(allMutEffects, env == chosenEnv)

# ~~ Step 2: Calculate regressions ####
# Choose variable slopes model since will just focus on results from 1 env
df_mutSlopes_concatData = allMutEffects %>%
  nest(data = -c(Mut_ID)) %>%
  mutate(fit = map(data, ~lm(avgS~Mean_GR, .)), glance = map(fit, ~tidy(.))) %>%
  unnest(glance) %>%
  dplyr::select( Mut_ID, term,estimate, std.error, statistic, p.value ) %>%
  pivot_wider(names_from = term, 
              values_from = c(estimate, std.error, statistic, p.value),
              names_sep = "_") %>%
  dplyr::select( Mut_ID,`estimate_(Intercept)` ,estimate_Mean_GR, `p.value_(Intercept)`, p.value_Mean_GR)

names(df_mutSlopes_concatData) = c( 'Mut_ID', 'intercept', 'slope', 'p_val_int', 'p_val_slope')




# ~~~ Plot 1: Hist of concat data slopes ####
# E.g., slopes of growth rate versus fitness effect 
# for data including both tn-insertions (directly measured)
# and reversions (infered)

dist_concatDataSlopes = ggplot(df_mutSlopes_concatData, aes(x = slope))+
  stat_binline(binwidth = 0.1, aes(y = after_stat(density)*0.1),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_x_continuous(limit = c(-1.1,0.1))+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.5, lineend = 'round')+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  #ggtitle('All Regressions')+
  xlab('Slope')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')

#ggsave(paste0(figSave, 'dist_concatDataSlopes.pdf'),dist_concatDataSlopes, unit = 'cm', width =4, height = 6 )



sum(df_mutSlopes_concatData$slope<0)/length(df_mutSlopes_concatData$slope)


### TAKEAWAY: 100% of slopes are negative in the concatonated data set 



# ~~~ Plot 2: Plot all regressions concat data ####
allReg_ConcatData = ggplot()+
  geom_abline(data = subset(df_mutSlopes_concatData), aes(slope = slope, intercept = intercept, group = Mut_ID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.9,lineend='round')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
  geom_vline(xintercept = grShift, color = 'grey', linewidth = 1, linetype = 'solid')+
  scale_color_manual( values = c('black', 'black', 'black'))+
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(-0.4,0.2))+
  xlab('Growth Rate (1/hr)')+
  ylab('Fitness effect')+
  geom_text(aes(label = 'Data', x = 0.05, y = -0.4), color = 'black', size = 3)+
  #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.title = element_text(color = 'black', size = 8),
        axis.text =  element_text(color = 'black', size = 8),
        legend.text = element_text(color = 'black', size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave, 'allReg_ConcatData.pdf'),allReg_ConcatData, unit = 'cm',width = 6.3, height = 3)



### TAKEAWAY: See clear 'bow-tie' plot among the concatenated data which includes both insertions and reversions



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Parameterize genotype to fitness map ####
# We now calculate mu_i and v_i values, as described in our 
# supplementary materials and in Reddy_desai eLife 2021
# The main takeaway here is that the only input into the model 
# is the distribution of global epistasis slopes from the concatonated data



### ~~~~~




## Create  4 site landscape
allSlopes_data = df_mutSlopes_concatData$slope

numSites = 4
idsChoice = seq(1,length(allSlopes_data), length.out = numSites+1)
idsChoice = floor(idsChoice)[1:numSites]
slopesChoice = allSlopes_data[order(allSlopes_data)]
slopesChoice = slopesChoice[idsChoice]

# Parameterize model ####
vis = -1/2 * slopesChoice
mus = vis/ ( 1 - vis)

mu = mean(mus) # mean mu of DVF (see reddy and desai)
p = length(mus) # number of loci
M = 200 # number of independent pathways -> they use 400 loci and 500 paths. I have only 94 paths, so use 200 paths. More paths = more computationally intensive 
maxn_m = numSites # max epistatic order = numsites in small landscape
modules = generate_modules(mus, M, maxn_m) # creates weights and loci
weights = modules[[1]]    # weights is list of all the epistasis coefficients (sampled from standard normal) for each module 
loci = modules[[2]]   # loci is list of the loci IDs in each module (each has between 2 and max_m loci)
Ks_all = get_Ks_all(maxn_m) # store a list for each possible set of loci IDS for any given number of interacting pathways 0-max epistatic order


#  ~~ Normalize variance ####

# Calculate n_ms
n_ms <- sapply(loci, length)

# Normalize variance
Vloci <- numeric(p)
V <- 0
for (m in 1:length(weights)) {
  n_m <- length(loci[[m]])
  Ks <- Ks_all[[n_m]] 
  V <- V + sum(weights[[m]][2:length(weights[[m]])]^2) 
  
  for (i in 1:n_m) {
    filt <- which(Ks[, i] == 1) # get all loci where valye is 1
    Vloci[loci[[m]][i]] <- Vloci[loci[[m]][i]] + sum(weights[[m]][filt]^2)
  }
}
Vloci <- Vloci / V

# Normalize weights
for (m in 1:length(weights)) {
  weights[[m]] <- weights[[m]] / sqrt(V)
}


## ~~ !! NOTE:: Directly ABOVE is a RANDOM processs which will create a new landscape each time it is run ####
# We show one particular landscape generated with this model in figure s13

# Get fitness of ALL genotypes ####

allGenotypes = expand.grid(rep(list(c(-1,1)), numSites))
allGenotypes$num1s = rowSums(allGenotypes == 1)

# add other axis counter
allGenotypes = allGenotypes %>%
  group_by(num1s) %>%
  mutate(counter_this_1num = 1:length(num1s))
allGenotypes$counter_this_1num[allGenotypes$num1s %in% c(0,4)] =  allGenotypes$counter_this_1num[allGenotypes$num1s %in% c(0,4)]+2.5
allGenotypes$counter_this_1num[allGenotypes$num1s %in% c(1,3)] =  allGenotypes$counter_this_1num[allGenotypes$num1s %in% c(1,3)]+1


fit_backgrounds_allallGenotypes = get_fitness(weights,loci,allGenotypes, Ks_all)*grScale + grShift 

allGenotypes$Fitness = fit_backgrounds_allallGenotypes

allGenotypes$name = paste0(allGenotypes$Var1, allGenotypes$Var2,allGenotypes$Var3,allGenotypes$Var4)

fitLandscape = 
  ggplot(allGenotypes, aes(y = num1s, x = counter_this_1num))+
  geom_rect(aes(ymin = num1s-0.3, xmin =counter_this_1num-0.46,
                ymax = num1s+0.09, xmax = counter_this_1num + 0.46,
                fill = Fitness),
            color = 'black', linewidth = 1.1 )+
  geom_rect(aes(ymin = num1s-0.07, xmin =counter_this_1num-0.42,
                ymax = num1s+0.07, xmax = counter_this_1num + 0.42,
                fill = Fitness),
            color = 'white', fill = 'white')+
  #scale_fill_viridis()+
  scale_fill_gradient2(low = 'white', high = 'black')+
  # geom_text(size = 4, aes(label = name), color = 'black')+
  geom_text(size = 4, aes(label = Var1, x = counter_this_1num-0.22), color = '#B92A20')+
  geom_text(size = 4, aes(label = Var2, x = counter_this_1num-0.073), color = '#36CAE0')+
  geom_text(size = 4, aes(label = Var3, x = counter_this_1num+0.073), color = '#F244BB')+
  geom_text(size = 4, aes(label = Var4, x = counter_this_1num+0.22), color = '#FFC11F')+
  geom_text(size = 5, aes(label = round(Fitness,3), y = num1s-0.2), color = 'white')+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(size = 20, color = 'black'),
        axis.text =element_text(size = 18, color = 'black', family = 'Helvetica'),
        legend.position = 'top',
        plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
        axis.ticks.length=unit(.05, "cm"))


#ggsave(paste0('~/Dropbox/SKLab/Projects/TnSeq/Revision_1_ALL/_forOlivier/fitLandscape_greyscale.pdf'),fitLandscape, width = 20, height = 20, unit = 'cm')


## Global epi
sitechoice =1

allplots = list()


for(sitechoice in 1:numSites)
{
  allbgs_fit  = get_fitness(weights,loci,allGenotypes, Ks_all)*grScale + grShift 
  
  allmuts = allGenotypes
  allmuts[,sitechoice]  = -1*allmuts[,sitechoice]
  allmuts_fit  = get_fitness(weights,loci,allmuts, Ks_all)*grScale + grShift 
  
  alldeltaF = allmuts_fit-allbgs_fit
  
  dftest = data.frame(allbgs_fit, alldeltaF)
  dftest$orig1 = allGenotypes[,sitechoice] == 1
  
  dftest$orig1[ dftest$orig1 == T] = 'geneLoss' # goes from 1 -> -1 
  dftest$orig1[ dftest$orig1 == F] = 'geneGain' # goes from -1 -> 1
  
  
  
  allplots[[sitechoice]] = ggplot(dftest, aes(x = allbgs_fit,y = alldeltaF, group = orig1, color = orig1 ))+
    geom_point(size = 0.5)+
    scale_color_manual(values = c('orange', 'blue'))+
    geom_smooth(method = 'lm', se = F, linewidth = 0.5, lineend = 'round')+
    xlab('Background strain Fitness, 1/h')+
    ylab('Fitness effect, 1/h')+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.3, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.3, lineend = 'round'), 
          axis.title = element_blank(),
          axis.text =element_text(size = 10, color = 'black', family = 'Helvetica'),
          legend.position = 'right',
          plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
          axis.ticks.length=unit(.05, "cm"))
}


allGlobalEpi = ggarrange(plotlist = allplots, nrow = 2, ncol = 2, common.legend = T)
#ggsave(paste0('~/Dropbox/SKLab/Projects/TnSeq/Revision_1_ALL/_forOlivier/allGlobalEpi_XYZ.pdf'),allGlobalEpi, width = 9, height = 11, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S13: Full Reddy Desai Implementation in R ####


# SETUP ####

# ~~ Set GR shift and Scale ####
# ~~~ NOTE: these values define how I convert 'fitness' values from the model to realistic fitness values
grScale = 0.2 # usually 0.2 to produce Rsq of slope ~ int = 0.98
grShift = 0.245 # equivalent to the pivot GR, same as 30sc5 in concatenated data, found by running code below then adjusting accordingly
chosenEnv = '30SC5' # choose one environment to do concatenated data regressions with 







# ~~~ Plot 1: Hist of concat data slopes ####
# E.g., slopes of growth rate versus fitness effect 
# for data including both tn-insertions (directly measured)
# and reversions (infered)

dist_concatDataSlopes = ggplot(df_mutSlopes_concatData, aes(x = slope))+
  stat_binline(binwidth = 0.1, aes(y = after_stat(density)*0.1),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_x_continuous(limit = c(-1.1,0.1))+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.5, lineend = 'round')+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  #ggtitle('All Regressions')+
  xlab('Slope')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')

#ggsave(paste0(figSave, 'dist_concatDataSlopes.pdf'),dist_concatDataSlopes, unit = 'cm', width =4, height = 6 )



sum(df_mutSlopes_concatData$slope<0)/length(df_mutSlopes_concatData$slope)


### TAKEAWAY: 100% of slopes are negative in the concatonated data set 



# ~~~ Plot 2: Plot all regressions concat data ####
allReg_ConcatData = ggplot()+
  geom_abline(data = subset(df_mutSlopes_concatData), aes(slope = slope, intercept = intercept, group = Mut_ID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.9,lineend='round')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
  geom_vline(xintercept = grShift, color = 'grey', linewidth = 1, linetype = 'solid')+
  scale_color_manual( values = c('black', 'black', 'black'))+
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(-0.4,0.2))+
  xlab('Growth Rate (1/hr)')+
  ylab('Fitness effect')+
  geom_text(aes(label = 'Data', x = 0.05, y = -0.4), color = 'black', size = 3)+
  #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.title = element_text(color = 'black', size = 8),
        axis.text =  element_text(color = 'black', size = 8),
        legend.text = element_text(color = 'black', size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave, 'allReg_ConcatData.pdf'),allReg_ConcatData, unit = 'cm',width = 6.3, height = 3)



### TAKEAWAY: See clear 'bow-tie' plot among the concatenated data which includes both insertions and reversions



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Parameterize genotype to fitness map ####
# We now calculate mu_i and v_i values, as described in our 
# supplementary materials and in Reddy_desai eLife 2021
# The main takeaway here is that the only input into the model 
# is the distribution of global epistasis slopes from the concatonated data


df_mutSlopes_concatData$vi = -1/2 * df_mutSlopes_concatData$slope 
df_mutSlopes_concatData$mu_i = df_mutSlopes_concatData$vi / ( 1 - df_mutSlopes_concatData$vi) # note this requires the data to only have 1 env, since in this case it only has1 row per site
df_mutSlopes_concatData$mutID_reCount = 1:nrow(df_mutSlopes_concatData) # include new mutID because model has just the mus in order, so lose the actual Mut_ID mapping

mus = df_mutSlopes_concatData$mu_i
mu = mean(mus) # mean mu of DVF (see reddy and desai)
p = length(mus) # number of loci
M = 200 # number of independent pathways -> they use 400 loci and 500 paths. I have only 94 paths, so use 200 paths. More paths = more computationally intensive 
maxn_m = 14 # max epistatic order
modules = generate_modules(mus, M, maxn_m) # creates weights and loci
weights = modules[[1]]    # weights is list of all the epistasis coefficients (sampled from standard normal) for each module 
loci = modules[[2]]   # loci is list of the loci IDs in each module (each has between 2 and max_m loci)
Ks_all = get_Ks_all(maxn_m) # store a list for each possible set of loci IDS for any given number of interacting pathways 0-max epistatic order
# used to efficiently determine whether a given epistatic weight will have a + or - contribution for a given locus
# when implementing eq(1) from Reddy and Desai 




## ((IMPORTANT)): Sampling the modules is a random process. This section above will define a 
# SINGLE, well-behaved genotype to fitness map. The rest of the script
# explores global epistasis and evolutionary outcomes on this map
# We report one instance of a map made with our model in S13


#  ~~ Normalize variance ####

# Calculate n_ms
n_ms <- sapply(loci, length)

# Normalize variance
Vloci <- numeric(p)
V <- 0
for (m in 1:length(weights)) {
  n_m <- length(loci[[m]])
  Ks <- Ks_all[[n_m]] 
  V <- V + sum(weights[[m]][2:length(weights[[m]])]^2) 
  
  for (i in 1:n_m) {
    filt <- which(Ks[, i] == 1) # get all loci where valye is 1
    Vloci[loci[[m]][i]] <- Vloci[loci[[m]][i]] + sum(weights[[m]][filt]^2)
  }
}
Vloci <- Vloci / V

# Normalize weights
for (m in 1:length(weights)) {
  weights[[m]] <- weights[[m]] / sqrt(V)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Get distribution of random strain GRs on this G-P landscape####


# First, sample many random individuals and calculate fitness
# (WARNING) takes a couple of minutes 
numBgs_i = 5000 # number backgrounds 
large_randomStrainSet =   2 * (matrix(rnorm(numBgs_i * p), numBgs_i, p) > 0) - 1 # sample many individuals with 'genotype' = sequence of 1/-1 values
fit_backgrounds_largeSet = get_fitness(weights,loci,large_randomStrainSet, Ks_all)*grScale + grShift # calculate unadjusted fitness of all individuals 
order_fits_bgs = order(fit_backgrounds_largeSet)


# ~~~ Plot 3: Plot hist of large GR set ####
dffits = data.frame(fit = fit_backgrounds_largeSet)

# Get fitness of the all [1] genotype and the all [-1] genotype
oneIndiv = t(matrix(numeric(100)))+1 
fit1= get_fitness(weights,loci,oneIndiv, Ks_all)*grScale + grShift


#write.csv(dffits, paste0(fileSave, 'bgStrainGRs.csv'))

histStrainGrs = ggplot(dffits, aes(x = fit))+
  stat_binline(binwidth = 0.1, aes(y = after_stat(density)*0.1),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  geom_histogram(binwidth = 0.1, aes(y = after_stat(density)*0.1), fill = 'grey' , color = 'grey')+
  scale_x_continuous(limit = c(-0.5,1.2))+
  geom_vline(xintercept = grShift, linetype = 'solid', color = 'white', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = fit1, linetype = 'solid', color = 'palegreen4', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = 0.02, linetype = 'solid', color = 'purple4', linewidth = 0.75, lineend = 'round')+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  xlab('Fitness')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')

#ggsave(paste0(figSave,'histStrainGrs_RD.pdf'), histStrainGrs, width = 6, height = 7, unit = 'cm')



# TAKEAWAY: The distribution of random strain fitnesses is a normal distribution centered at the pivot GR (gr Shift value)
# The all [1] is neither the best nor the least fit







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~ Correlation Decay ####

numBgs_small = 100 # number backgrounds 
randomStrainSet_small =   2 * (matrix(rnorm(numBgs_small * p), numBgs_small, p) > 0) - 1 # sample many individuals with 'genotype' = sequence of 1/-1 values
fit_backgrounds_small = get_fitness(weights,loci,randomStrainSet_small, Ks_all)*grScale + grShift # calculate unadjusted fitness of all individuals 


df_smallerSet = data.frame(randomStrainSet_small)
df_smallerSet$name <- apply(df_smallerSet, 1, function(row) paste(row, collapse = ""))

df_smallerSet$Fitness = fit_backgrounds_small


library(stringdist)

# replace -1 with 0 
df_smallerSet$name_modified <- gsub("-1", "0", df_smallerSet$name)
df_smallerSet$ID = 1



sample_vectors <- function(original_vector, num_samples=30, sample_length=20) {
  samples <- matrix(nrow=num_samples, ncol=sample_length)
  for (i in 1:num_samples) {
    samples[i, ] <- sample(original_vector, size=sample_length, replace=FALSE)
  }
  return(samples)
}



# ~~ pick 3 starting bgs (low, high, near pivot)
minGR_bg = randomStrainSet_small[which.min(abs(fit_backgrounds_small-0.02)), ] # pick low positive GR
midGR_bg = randomStrainSet_small[which.min(abs(fit_backgrounds_small-grShift)), ]
maxGR_bg = randomStrainSet_small[which.max(fit_backgrounds_small), ]

startBgList = list(minGR_bg , midGR_bg, maxGR_bg)
names_bg = c('min', 'mid', 'max')

allNeighborSummary = data.frame()
for(focalBg_it in 1:3)
{
  print(paste0('BG: ', focalBg_it))
  focalBg = startBgList[[focalBg_it]]
  allOutDf = list()
  
  for(locusChoice in 1:94)
  {
    print(paste0('locus: ', locusChoice))
    
    
    
    ## get initial effects
    focalGR_bg_mut = focalBg
    focalGR_bg_mut[locusChoice] = -1*focalGR_bg_mut[locusChoice] 
    initialFitEffect =  (get_fitness(weights, loci, t(focalGR_bg_mut), Ks_all  )*grScale + grShift) -   (get_fitness(weights, loci, t(focalBg), Ks_all  )*grScale + grShift)
    
    
    
    # ~~ For each bg, generate many single mut neighbors, excluding the focal locus 
    numMutsToDo = 93
    numMutsTot = p
    
    allDfout = data.frame()
    for(distanceUse in 1:20)
    {
      print(distanceUse)
      mutSites = sample_vectors((1:numMutsTot)[-locusChoice],  num_samples=numMutsToDo, sample_length=distanceUse )
      
      # rep focal bg before mut,, make sure re-establish before each distance change 
      allBgs_Here = as.data.frame(matrix(rep(focalBg, times = numMutsToDo), ncol = numMutsTot, byrow = TRUE))
      
      
      
      # get mutants 
      for(i in 1:numMutsToDo) 
      {
        sites = mutSites[i, ]
        allBgs_Here[i, sites] = -1*allBgs_Here[i, sites] 
      }
      
      Bg_wFocalMut = allBgs_Here
      Bg_wFocalMut[, locusChoice ] = -1* Bg_wFocalMut[i, locusChoice] 
      
      
      
      # get fitness effect of mut in these bgs
      fitnessMuts = get_fitness(weights, loci, allBgs_Here, Ks_all  )*grScale + grShift
      fitnessMuts_wFocal = get_fitness(weights, loci, Bg_wFocalMut, Ks_all  )*grScale + grShift
      
      fitEffect = fitnessMuts_wFocal - fitnessMuts
      
      
      df_out = data.frame(origFitEffect = initialFitEffect, newFitEffect = fitEffect, dist = distanceUse, locusID = locusChoice)
      allDfout = rbind(allDfout,df_out )
      
      
    }
    
    allOutDf[[locusChoice]] = allDfout
    
  }
  
  allOutDf_unlist =  bind_rows(allOutDf, .id = "ID")
  
  
  summary_neighborhoodCor = allOutDf_unlist %>%
    group_by(dist) %>%
    summarize(corEffects = cor(origFitEffect, newFitEffect))
  
  
  summary_neighborhoodCor$typeStart = names_bg[focalBg_it]
  allNeighborSummary = rbind(allNeighborSummary,summary_neighborhoodCor )
  
}


neighborCorDecay = ggplot(allNeighborSummary, aes(x = dist, y = corEffects, group = typeStart, color = typeStart))+
  geom_point(shape = 16)+
  geom_line(linewidth = 0.5)+
  xlab('Genotypic distance')+
  ylab("Mutation effect correlation")+
  scale_color_manual(values = c('#802922', '#00BFAE', '#BFBE9C'))+
  geom_hline(yintercept = 0, color ='grey')+
  scale_y_continuous(trans = 'log', breaks = c(0.05, 0.15,0.45,1))+
  geom_abline(slope = -0.15, intercept = 0, linewidth = 1.2)+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.position = 'right')

#ggsave(paste0('~/Dropbox/SKLab/Projects/TnSeq/Revision_1_ALL/Figures_Revision/neighborCorDecay.pdf'),neighborCorDecay, unit = 'cm', width = 7, height = 7)

allNeighborSummary$log10_y = log(allNeighborSummary$corEffects)


summary(lm(log10_y~dist,data = subset(allNeighborSummary, typeStart == 'min')))
summary(lm(log10_y~dist,data = subset(allNeighborSummary, typeStart == 'mid')))
summary(lm(log10_y~dist,data = subset(allNeighborSummary, typeStart == 'max')))




## add hist of fits

testFit_1 = get_fitness(weights,loci,t(minGR_bg), Ks_all)*grScale + grShift
testFit_2 = get_fitness(weights,loci,t(midGR_bg), Ks_all)*grScale + grShift
testFit_3 = get_fitness(weights,loci,t(maxGR_bg), Ks_all)*grScale + grShift

#write.csv(dffits, paste0(fileSave, 'bgStrainGRs.csv'))

histStrainGrs2 = ggplot(dffits, aes(x = fit))+
  stat_binline(binwidth = 0.1, aes(y = after_stat(density)*0.1),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  geom_histogram(binwidth = 0.1, aes(y = after_stat(density)*0.1), fill = 'grey' , color = 'grey')+
  scale_x_continuous(limit = c(-0.5,1.2))+
  geom_vline(xintercept = grShift, linetype = 'solid', color = 'white', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = fit1, linetype = 'solid', color = 'palegreen4', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = 0.02, linetype = 'solid', color = 'purple4', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = testFit_1, linetype = 'solid', color = '#BFBE9C', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = testFit_2, linetype = 'solid', color = '#00BFAE', linewidth = 0.75, lineend = 'round')+
  geom_vline(xintercept = testFit_3, linetype = 'solid', color = '#802922', linewidth = 0.75, lineend = 'round')+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  xlab('Fitness')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 10),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')

#ggsave(paste0('~/Dropbox/SKLab/Projects/TnSeq/Revision_1_ALL/Figures_Revision/histStrainGrs_wCorDecayPts.pdf'),histStrainGrs2,width = 6, height = 7, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Get fitness effects at all loci across many bgs ####

# sample down to 300 backgrounds evenly spaced in fitness
numBgs = 300 ## Actual number to use in global epistasis check
idsTopick = ceiling(seq(1,numBgs_i, length.out = numBgs)) # round to integer
large_randomStrainSet = large_randomStrainSet[order_fits_bgs, ]
smallerStrainSet = large_randomStrainSet[idsTopick, ]
fit_backgrounds = get_fitness(weights,loci,smallerStrainSet, Ks_all)*grScale + grShift

# Sanity check: re-calculated fitnesses are same as already had
# plot(fit_backgrounds_largeSet[order_fits_bgs][idsTopick], fit_backgrounds)


## Calculate 'fitness effect' of flipping each mutation (e.g., -1 -> 1 or 1-> -1)
# at each site in all backgrounds 
# (WARNING) takes a few minutes, counts the locus it is on (goes to 94)
df_s_ests_RD_list = list()
for( pickMut in 1:p)
{
  print(pickMut)
  slopeLoci =  df_mutSlopes_concatData$slope[pickMut]
  intLoci =  df_mutSlopes_concatData$intercept[pickMut]
  
  
  mutBgs = smallerStrainSet
  mutBgs[, pickMut] = -1*smallerStrainSet[,pickMut]
  
  fit_muts = get_fitness(weights,loci,mutBgs, Ks_all)*grScale + grShift
  
  df_mut = data.frame(fitness_bg_RD = fit_backgrounds, fitness_mut_RD = fit_muts)
  df_mut$mutID =  df_mutSlopes_concatData$mutID_reCount[pickMut]
  df_mut$predSlope = slopeLoci
  df_mut$predInt = intLoci
  
  
  df_s_ests_RD_list[[pickMut]] = df_mut
  # ggplot(df_mut, aes(x = fitness_bg_RD, y = fitness_mut_RD -fitness_bg_RD ))+
  #   geom_point(shape = 16)+
  #   geom_abline(slope = slopeLoci, intercept = 0, col = 'red')+
  #   # ylim(-0.1,0.02)+
  #   geom_smooth(method = 'lm', se = F)+
  #   theme_classic()
}

df_s_ests_RD = bind_rows(df_s_ests_RD_list)
df_s_ests_RD$fitEffect = df_s_ests_RD$fitness_mut_RD - df_s_ests_RD$fitness_bg_RD

#write.csv(df_s_ests_RD, paste0(fileSave,'df_s_ests_RD.csv'))




## CHECK that the global epistasis slopes and intercepts we calculate from the genotype-fitness map we have
## match with the inputted desired slopes 

## ~~ Plot 3-4: Predicted v Actual Slope and Intercept ####
getslopes = df_s_ests_RD %>%
  nest(data = -c(mutID,predSlope,predInt)) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(fitEffect ~ fitness_bg_RD , .)), glance = map(fit, ~ tidy(.x))) %>%
  unnest(glance)


toSave = getslopes %>%
  dplyr::select(mutID, predSlope, predInt, term, estimate, std.error, statistic, p.value)
#write.csv(toSave, paste0(fileSave,'allSlopes_pred_actual_RD.csv'))


test = subset(getslopes, term == 'fitness_bg_RD') # get slope
test_slopesOnly = test[,c('mutID','predSlope', 'estimate')]
names(test_slopesOnly) = c('mutID','predSlope', 'slope')


test = subset(getslopes, term == '(Intercept)') # get intercept
test_intsOnly = test[,c('mutID','predInt', 'estimate')]
names(test_intsOnly) = c('mutID','predInt', 'intercept')

## ~~ see if slopes found here relate to the slopes assigned for each locus

pred_v_actualSlope = ggplot(test_slopesOnly, aes(x = predSlope, y = slope))+
  geom_point(shape = 16, alpha = 0.6)+
  geom_abline(slope = 1, intercept = 0, color ='grey', lineend = 'round', linewidth = 0.5)+
  xlim(-0.9,0.1)+
  ylim(-0.9,0.1)+
  ylab('Slope in G-F map')+
  xlab('Input slope')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        #axis.ticks.length=unit(.04, "cm"),
        axis.title = element_text(size=10,color = 'black'),
        axis.text = element_text(size=8,color = 'black'),
        legend.position = 'right' ,
        plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) 



pred_v_actualInt = ggplot(test_intsOnly, aes(x = predInt, y = intercept))+
  geom_point(shape = 16, alpha = 0.6)+
  geom_abline(slope = 1, intercept = 0, color ='grey', lineend = 'round', linewidth = 0.5)+
  xlim(-0.05,0.2)+
  ylim(-0.05,0.2)+
  ylab('Intercept in G-F map')+
  xlab('Estimated intercept')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        #axis.ticks.length=unit(.04, "cm"),
        axis.title = element_text(size=10,color = 'black'),
        axis.text = element_text(size=8,color = 'black'),
        legend.position = 'right' ,
        plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) 


summary(lm(slope~predSlope,data = test_slopesOnly))
summary(lm(intercept~predInt,data = test_intsOnly))

# ggsave(paste0(figSave, 'pred_v_actualSlope_RD.pdf'), pred_v_actualSlope, unit = 'cm', width = 3.3, height = 6)
# ggsave(paste0(figSave, 'pred_v_actualInt_RD.pdf'), pred_v_actualInt, unit = 'cm', width = 3.3, height = 6)


## TAKEAWAY: Slopes and intercepts in the G-F landscape are very similar to those
# we see in our data. Reminder, the only input parameter into the landscape was the slope
# of the microscopic global epistasis trends



# ~~ Test: Slope Intercept Correlation 
slopeInt_RD = getslopes %>%
  dplyr::select(mutID, predSlope, predInt, term, estimate)%>%
  pivot_wider(names_from = term, 
              values_from = c(estimate),
              names_sep = "_") %>%
  dplyr::select( mutID, predSlope, predInt,`(Intercept)` ,fitness_bg_RD)
names(slopeInt_RD)[4:5] = c('intercept', 'slope')


slopeIntCorr_RD = ggplot(slopeInt_RD, aes(x = slope, y = intercept))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_point( alpha = 0.5, size = 2, shape = 16)+
  geom_smooth(method = 'lm',formula = 'y~0+x', se = F, lineend = 'round')+
  xlab('Slope') +  ylab('Intercept')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 12),
        axis.text =  element_text(color = 'black', size = 10),
        legend.text = element_text(color = 'black', size = 10),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)


#ggsave(paste0(figSave,'slopeIntCorr_RD.pdf' ),slopeIntCorr_RD, width = 5, height = 5.4, unit = 'cm')


## ~~ Plot 5: Model bowtie plot ####
allReg_RD = ggplot()+
  geom_abline(data = subset(slopeInt_RD), aes(slope = slope, intercept = intercept, group = mutID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.9,lineend='round')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
  geom_vline(xintercept = grShift, color = 'grey', linewidth = 0.5, linetype = 'solid')+
  scale_color_manual( values = c('black', 'black', 'black'))+
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(-0.4,0.2))+
  xlab('Growth Rate (1/hr)')+ 
  ylab('Fitness effect')+
  geom_text(aes(label = 'Model', x = 0.05, y = -0.4), color = 'black', size = 3)+
  #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black',lineend='round'),
        axis.title = element_text(color = 'black', size = 8),
        axis.text =  element_text(color = 'black', size = 8),
        legend.text = element_text(color = 'black', size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(figSave, 'allReg_RD.pdf'),allReg_RD, unit = 'cm',width = 6.3, height = 3)


## TAKEAWAY: The G-F map also creates a bowtie plot which looks very much like the one we see in the data 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Insertions vs reversions in model ####

# ~~~~~ A. Single Examples
numBgs_i = 500 # sample many strains to get a good range above 0 
random_1 =   2 * (matrix(rnorm(numBgs_i * p), numBgs_i, p) > 0) - 1 
random_2 =   2 * (matrix(rnorm(numBgs_i * p), numBgs_i, p) > 0) - 1 

fits1 = get_fitness(weights,loci,random_1, Ks_all)*grScale + grShift
fits2 = get_fitness(weights,loci,random_2, Ks_all)*grScale + grShift

# reorder
random_1 = random_1[order(fits1), ]
fits1 = fits1[order(fits1) ]
random_2 = random_2[order(fits2), ]
fits2 = fits2[order(fits2) ]

# get rid of negative GRs
random_1 = random_1[fits1>0,]
fits1 = fits1[fits1>0]
random_2 = random_2[fits2>0,]
fits2 = fits2[fits2>0]


# pick 42 random along traj
toPick_1 = floor(seq(1, length(fits1), length.out = 42))
toPick_2 = floor(seq(1, length(fits2), length.out = 42))

random_1 = random_1[toPick_1, ]
random_2 = random_2[toPick_2, ]


# Get fitness effects for all of them at the focal locus

mutPick = 72 
# ~~ insertion (1 -> -1)
random_1_w1 = random_1
random_1_w1[, mutPick] = 1

bgFits = get_fitness(weights,loci,random_1_w1, Ks_all)*grScale + grShift

random_1_wMinus1 = random_1_w1
random_1_wMinus1[, mutPick] = -1
mutFits = get_fitness(weights,loci,random_1_wMinus1, Ks_all)*grScale + grShift
fitEffects = mutFits - bgFits
df_insertion = data.frame(bgFits, mutFits, fitEffects, type = 'insertion')


# ~~ reversion (-1 -> 1)
random_2_w1 = random_2
random_2_w1[, mutPick] = -1

bgFits = get_fitness(weights,loci,random_2_w1, Ks_all)*grScale + grShift
random_2_wMinus1 = random_2_w1
random_2_wMinus1[, mutPick] = 1
mutFits = get_fitness(weights,loci,random_2_wMinus1, Ks_all)*grScale + grShift
fitEffects = mutFits - bgFits
df_reversion = data.frame(bgFits, mutFits, fitEffects, type = 'reversion')


df_RD_allType_ex = rbind(df_insertion,df_reversion )
#write.csv(df_RD_allType_ex, paste0(fileSave, 'df_RD_allType_ex.csv'))
df_RD_allType_ex$type = as.factor(df_RD_allType_ex$type)


p1_insertion_reversion_RD = ggplot(df_RD_allType_ex, aes(x = bgFits, y = fitEffects, color = type))+
  geom_hline(yintercept = 0, color = 'grey', lineend = 'round', linewidth = 0.5)+
  geom_point(shape = 16)+
  geom_smooth(method = 'lm',se = F,linetype = 'solid', lineend = 'round', linewidth = 0.5)+
  ylab('Proportion')+
  scale_color_manual(values = c('#56239E', '#E4B02D'), drop = F)+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        #axis.ticks.length=unit(.04, "cm"),
        axis.title = element_text(size=10,color = 'black'),
        axis.text = element_text(size=8,color = 'black'),
        legend.position = 'right' ,
        plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) 




# ~~~~~ B. Distribution of slopes
## ~~ Just use same random genotypes from above
allDF_RD_allType_list = list()
allSlope_int_allType_RD_list = list()
for(mutPick in 1:p)
{
  print(mutPick)
  # ~~ insertion (1 -> -1)
  random_1_w1 = random_1
  random_1_w1[, mutPick] = 1
  
  bgFits = get_fitness(weights,loci,random_1_w1, Ks_all)*grScale + grShift
  
  random_1_wMinus1 = random_1_w1
  random_1_wMinus1[, mutPick] = -1
  mutFits = get_fitness(weights,loci,random_1_wMinus1, Ks_all)*grScale + grShift
  fitEffects = mutFits - bgFits
  df_insertion = data.frame(bgFits, mutFits, fitEffects, type = 'insertion')
  
  # ~~ reversion (-1 -> 1)
  random_2_w1 = random_2
  random_2_w1[, mutPick] = -1
  
  bgFits = get_fitness(weights,loci,random_2_w1, Ks_all)*grScale + grShift
  random_2_wMinus1 = random_2_w1
  random_2_wMinus1[, mutPick] = 1
  mutFits = get_fitness(weights,loci,random_2_wMinus1, Ks_all)*grScale + grShift
  fitEffects = mutFits - bgFits
  df_reversion = data.frame(bgFits, mutFits, fitEffects, type = 'reversion')
  
  
  df_RD_allType_thismut = rbind(df_insertion,df_reversion )
  #write.csv(df_RD_allType_thismut, paste0(fileSave, 'df_RD_allType_thismut_', mutPick, '.csv'))
  df_RD_allType_thismut$type = as.factor(df_RD_allType_thismut$type)
  df_RD_allType_thismut$mutID = mutPick
  
  
  
  ## Get slope and int
  slopesInts_mutTypes_RD = df_RD_allType_thismut %>%
    nest(data = -type) %>%
    mutate(fit = map(data, ~lm(fitEffects~bgFits, data = .)), tidy= map(fit, ~tidy(.))) %>%
    unnest(tidy) %>%
    dplyr::select(type, term, estimate) %>%
    pivot_wider(names_from = term, 
                values_from = estimate)
  
  names(slopesInts_mutTypes_RD) = c('type', 'intercept', 'slope')
  slopesInts_mutTypes_RD$mutID = mutPick
  
  allDF_RD_allType_list[[mutPick]] = df_RD_allType_thismut
  allSlope_int_allType_RD_list[[mutPick]] = slopesInts_mutTypes_RD
  
}

allDF_RD_allType = bind_rows(allDF_RD_allType_list)
allSlope_int_allType_RD = bind_rows(allSlope_int_allType_RD_list )
allSlope_int_allType_RD$type = as.factor(allSlope_int_allType_RD$type)


p3_Hist_insertion_reversion_RD = ggplot(allSlope_int_allType_RD, aes(x = slope, color = type))+
  stat_binline(binwidth = 0.05, aes(y = after_stat(density)*0.05),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  ylab('Proportion')+
  scale_x_continuous(limits = c(-1,0.3))+
  scale_color_manual(values = c('#56239E', '#E4B02D'), drop = F)+
  geom_vline(xintercept = 0, color = 'grey')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        #axis.ticks.length=unit(.04, "cm"),
        axis.title = element_text(size=10,color = 'black'),
        axis.text = element_text(size=8,color = 'black'),
        legend.position = 'right' ,
        plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) 





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Evolution Simulations ####


# ~~ Part 1: Starting from [all 1] type ####
randIndiv =t(matrix(numeric(p)))+1 
print(get_fitness(weights, loci, randIndiv, Ks_all)*grScale+grShift) 


## (WARNING) takes ~2 min/iteration!! 
## for the text, we run 50 iterations of 50k gens each. 
## for demonstration purposes here, I just run 15 iterations of 5k gens each
nIts = 15
testSim = evolve_ngens(randIndiv, 5000, nIts, U = 10^-5, N = 10^4, weights, loci, Ks_all)


df_allIt = testSim[[1]]
# write.csv(df_allIt, paste0(fileSave, 'df_allIt_all1start.csv'))


df_allIt_avg = df_allIt %>%
  group_by( time) %>%
  summarize(meanFit = mean(fitness), meanDiffloci = mean(diffFromFounder),meanNum_neg1 =mean(numNeg1))


tfrep = nIts<10
chosen_exIts = sample( 1:nIts,10, replace = tfrep)
subdf_allIt = subset(df_allIt, it%in%chosen_exIts)

fitTraj_RD = ggplot()+
  geom_line(data = subdf_allIt, aes(x = time, y = fitness, group = it), alpha = 0.5, linewidth = 0.4)+
  geom_line(data =df_allIt_avg, aes(x = time, y = meanFit), alpha = 1, linewidth = 1, color ='palegreen4' )+
  theme_classic()+
  xlab('Time, in generations')+
  ylab('Fitness')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 8),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')


numSubTraj_RD = ggplot()+
  geom_line(data = subdf_allIt, aes(x = time, y = numNeg1, group = it), alpha = 0.5, linewidth = 0.4)+
  geom_line(data =df_allIt_avg, aes(x = time, y = meanNum_neg1), alpha = 1, linewidth = 1, color ='palegreen4' )+
  theme_classic()+
  xlab('Time, in generations')+
  ylab('# KOS')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 8),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')


# ggsave(paste0(figSave, 'fitTraj_RD_all1.pdf'), fitTraj_RD, unit = 'cm',  width = 5.7, height = 3.5)
# ggsave(paste0(figSave, 'numSubTraj_all1_RD.pdf'), numSubTraj_RD, unit = 'cm', width = 5.7, height = 3.5)





# ~~ Part 2: Starting from arbitrary low fit type ####

n = 1
randIndiv = 2 * (matrix(rnorm(n* p), n, p) > 0) - 1 # !!! run until get low positive (~0.03) fitness
print(get_fitness(weights, loci, randIndiv, Ks_all)*grScale+grShift) # 
randIndiv_lowGR = randIndiv

## make sure grShift and grScale are defines
nIts = 15
testSim2 = evolve_ngens(randIndiv_lowGR, 5000, nIts, U = 10^-5, N = 10^4, weights, loci, Ks_all)

df_allIt_lowFit = testSim2[[1]]
# write.csv(df_allIt, paste0(fileSave, 'df_allIt_lowFit.csv'))
# write.csv(randIndiv_lowGR, paste0(fileSave, 'randIndiv_lowGR.csv'))


df_allIt_avg_lowFit = df_allIt_lowFit %>%
  group_by( time) %>%
  summarize(meanFit = mean(fitness), meanDiffloci = mean(diffFromFounder),meanNum_neg1 =mean(numNeg1))


tfrep = nIts<10
chosen_exIts = sample( 1:nIts,10, replace = tfrep)
subdf_allIt_lowFit = subset(df_allIt_lowFit, it%in%chosen_exIts)

fitTraj_RD_low = ggplot()+
  geom_line(data = subdf_allIt_lowFit, aes(x = time, y = fitness, group = it), alpha = 0.5, linewidth = 0.4)+
  geom_line(data =df_allIt_avg_lowFit, aes(x = time, y = meanFit), alpha = 1, linewidth = 1, color ='purple4' )+
  theme_classic()+
  xlab('Time, in generations')+
  ylab('Fitness')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 8),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')



numSubTraj_RD_low = ggplot()+
  geom_line(data = subdf_allIt_lowFit, aes(x = time, y = numNeg1, group = it), alpha = 0.5, linewidth = 0.4)+
  geom_line(data =df_allIt_avg_lowFit, aes(x = time, y = meanNum_neg1), alpha = 1, linewidth = 1, color ='purple4' )+
  theme_classic()+
  xlab('Time, in generations')+
  ylab('# KOS')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 8),
        axis.text = element_text(color = 'black', size = 8),
        legend.position = 'right')

# ggsave(paste0(figSave, 'fitTraj_RD_low.pdf'), fitTraj_RD_low, unit = 'cm', width = 5.7, height = 3.5)
# ggsave(paste0(figSave, 'numSubTraj_RD_low.pdf'), numSubTraj_RD_low, unit = 'cm', width = 5.7, height = 3.5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####




