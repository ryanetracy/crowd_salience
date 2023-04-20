#############################
# crowding disease analyses
#############################

# list packages, check installation, and load
pckgs <- c(
  "tidyverse",
  "lme4",
  "lmerTest",
  "interactions",
  "effectsize",
  "performance",
  "haven",
  "effsize",
  "emmeans",
  "ggcorrplot",
  "EMAtools",
  "HLMdiag",
  "see"
  )

for (pckg in pckgs) {
  if (
    !(pckg %in% installed.packages())
    ) {
    install.packages(pckg)
  }
  lapply(pckg, library, character.only = T)
}

# data
cdData <- read_sav("Crowding Disease Face.sav")

# get a sense of the group sizes (crowd = 1, control = 2)
plyr::count(cdData$Condition)

# and the sources (MS = 1, NJ = 2)
plyr::count(cdData$Source)

# come up with a vector of new names 
newNames <- paste0(
  rep("stim"),
  rep(1:4, each = 2),
  rep(c("D1", "D2", "H1", "H2"), each = 2),
  rep(c("Hands", "Bus"))
)

# apply the new names to the data set
names(cdData)[28:35] <- newNames


# #  conduct analyses with 'affect' as a covariate
cdLongA <- cdData %>%
  select(
    id,
    Condition,
    Source,
    Affect,
    stim1D1Hands,
    stim1D1Bus,
    stim2D2Hands,
    stim2D2Bus,
    stim3H1Hands,
    stim3H1Bus,
    stim4H2Hands,
    stim4H2Bus
    ) %>%
  gather(key = "stimID", value = "response",
         stim1D1Hands:stim4H2Bus) %>%
  separate(col = "stimID", into = c("stimID", "outcome"), sep = 5) %>%
  separate(col = "outcome", into = c("health", "outcome"), sep = 2) %>%
  # make the 'hands' and 'bus' variables two different outcomes 
  pivot_wider(names_from = "outcome", values_from = "response") %>%
  # this is just to rename/recode some of the variables to make it easier to understand the contrasts later
  mutate(health = case_when(
    health == "D1" | health == "D2" ~ "Diseased",
    health == "H1" | health == "H2" ~ "Healthy"
  )) %>%
  mutate(Condition = case_when(
    Condition == 1 ~ "Crowd",
    Condition == 2 ~ "Control"
  )) %>%
  mutate(Source = case_when(
    Source == 1 ~ "MS",
    Source == 2 ~ "NJ"
  )) %>%
  # and now create linear contrast variables (highly recommended for lmms)
  mutate(condC = case_when(
    Condition == "Crowd" ~ 1,
    Condition == "Control" ~ -1
  )) %>%
  mutate(sourceC = case_when(
    Source == "MS" ~ 1,
    Source == "NJ"~ -1
  )) %>%
  mutate(healthC = case_when(
    health == "Diseased" ~ 1,
    health == "Healthy" ~ -1
  ))

# show the first 10 rows:
head(cdLongA, 10)


# start with boxplots
# for the 'hands' variable
ggplot(cdLongA, aes(health, Hands, fill = Condition)) +
  geom_boxplot(outlier.color = "black") +
  theme_classic(base_size = 15) +
  facet_grid(~ Source) +
  scale_fill_manual(values = c("gray70", "turquoise4"))

# and for the 'bus' variable
ggplot(cdLongA, aes(health, Bus, fill = Condition)) +
  geom_boxplot(outlier.color = "black") +
  theme_classic(base_size = 15) +
  facet_grid(~ Source) +
  scale_fill_manual(values = c("gray70", "turquoise4"))

# look at distributions (collapsed across conditions)
ggplot(cdLongA, aes(Hands)) +
  geom_histogram(bins = 7, fill = "turquoise4", alpha = .8, color = "black") +
  theme_classic(base_size = 15)

ggplot(cdLongA, aes(Bus)) +
  geom_histogram(bins = 7, fill = "turquoise4", alpha = .8, color = "black") +
  theme_classic(base_size = 15)



# run correlation analyses of the two outcome variables
# overall effects
cdParts <- cdLongA %>% 
  group_by(id) %>% 
  summarize(
    mHands = mean(Hands, na.rm = T),
    mBus = mean(Bus, na.rm = T)
    ) %>%
  as.data.frame
cor.test(cdParts$mHands, cdParts$mBus) # r(237) = .87, p < .001, r^2 = .76

# look at the items for each target category separately
cdPartsT <- cdLongA %>%
  group_by(id, health) %>%
  summarize(
    mHands = mean(Hands, na.rm = T),
    mBus = mean(Bus, na.rm = T)
    ) %>%
  as.data.frame
healthy <- filter(cdPartsT, health == "Healthy")
disease <- filter(cdPartsT, health == "Diseased")

cor.test(healthy$mHands, healthy$mBus)
cor.test(disease$mHands, disease$mBus)

# average the two DVs to make a 'contact' variable
cdLongA$contact <- rowMeans(cdLongA[,7:8], na.rm = T)


# explore the distributions for this averaged variable
ggplot(cdLongA, aes(health, contact, fill = Condition)) +
  geom_boxplot(outlier.color = "black") +
  theme_classic(base_size = 15) +
  facet_grid(~Source) +
  scale_fill_manual(values = c("gray70", "turquoise4"))

ggplot(cdLongA, aes(contact)) +
  geom_histogram(bins = 7, fill = "turquoise4", alpha = .8, color = "black") +
  theme_classic(base_size = 15)


# test different models
m1 <- lmer(contact ~ healthC * condC * sourceC + Affect + (1|id) + (1|stimID), data = cdLongA)
m2 <- lmer(contact ~ healthC * condC * sourceC + Affect + (healthC|id) + (1|stimID), data = cdLongA)
m3 <- lmer(contact ~ healthC * condC * sourceC + Affect + (1|id:Source) + (1|stimID), data = cdLongA)
m4 <- lmer(contact ~ healthC * condC * sourceC + Affect + (healthC|id:Source) + (1|stimID), data = cdLongA)
m5 <- lmer(contact ~ healthC * condC * sourceC + Affect + (healthC|id:Source) + (0 + healthC|stimID), data = cdLongA) # this one
m6 <- lmer(contact ~ healthC * condC * sourceC + Affect + (healthC|id) + (0 + healthC|stimID), data = cdLongA)

# compare models
anova(m1, m2, m3, m4, m5, m6)
compare_performance(m1, m2, m3, m4, m5, m6)
# go with m5, which accounts for the most random effects

summary(m5)
standardize_parameters(m5)
eta_squared(m5)
icc(m5)
r2(m5)
confint(m5)

# descriptives for each condition
# health
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = "health",
                 na.rm = T)

# condition
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = "Condition",
                 na.rm = T)

# source
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = "Source",
                 na.rm = T)


# explore two-way interactions: health x cond
# control condition 
modC2.1 <- lmer(contact ~ healthC + Affect + (healthC|id:Source) + (0 + healthC|stimID), data = filter(cdLongA, condC == -1))
summary(modC2.1)
standardize_parameters(modC2.1)
eta_squared(modC2.1)
icc(modC2.1)
r2(modC2.1)
confint(modC2.1)

# crowd condition
modC2.2 <- lmer(contact ~ healthC + Affect + (healthC|id:Source) + (0+healthC|stimID), data = filter(cdLongA, condC == 1))
summary(modC2.2)
standardize_parameters(modC2.2)
eta_squared(modC2.2)
icc(modC2.2)
r2(modC2.2)
confint(modC2.2)

# descriptives
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = c("health", "Condition"),
                 na.rm = T)


# health x source two-way interaction:
# NJ participants 
modC3.1 <- lmer(contact ~ healthC + Affect + (healthC|id:Source) + (0+healthC|stimID), data = filter(cdLongA, sourceC == -1))
summary(modC3.1)
standardize_parameters(modC3.1)
eta_squared(modC3.1)
icc(modC3.1)
r2(modC3.1)
confint(modC3.1)

# MS participants
modC3.2 <- lmer(contact ~ healthC + Affect + (healthC|id:Source) + (0+healthC|stimID), data = filter(cdLongA, sourceC == 1))
summary(modC3.2)
standardize_parameters(modC3.2)
eta_squared(modC3.2)
icc(modC3.2)
r2(modC3.2)
confint(modC3.2)

# descriptives
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = c("health", "Source"),
                 na.rm = T)


# condition x source interaction (moderate)
# NJ participants  
modC4.1 <- lmer(contact ~ condC + Affect + (1|id:Source) + (1|stimID), data = filter(cdLongA, sourceC == -1))
summary(modC4.1)
standardize_parameters(modC4.1)
eta_squared(modC4.1)
icc(modC4.1)
r2(modC4.1)
confint(modC4.1)

# crowd condition
modC4.2 <- lmer(contact ~ condC + Affect + (1|id:Source) + (1|stimID), data = filter(cdLongA, sourceC == 1))
summary(modC4.2)
standardize_parameters(modC4.2)
eta_squared(modC4.2)
icc(modC4.2)
r2(modC4.2)
confint(modC4.2)

# descriptives
Rmisc::summarySE(cdLongA,
                 measurevar = "contact",
                 groupvars = c("Source", "Condition"),
                 na.rm = T)


# make a summary of the data
cdSum <- Rmisc::summarySE(cdLongA,
                          measurevar = "contact",
                          groupvars = c("health", "Condition", "Source"),
                          na.rm = T)
cdSum

# rename facet labels
sourceLabs <- c("Northern Participants", "Southern Participants")
names(sourceLabs) <- c("NJ", "MS")

# plot it
# (with jittered values)
ggplot(cdSum, aes(health, contact, fill = Condition)) + 
  geom_point(data = cdLongA,
             aes(health, contact, color = Condition), 
             position = position_jitterdodge(.15, .15, .9),
             alpha = .5,
             size = .75) +
  geom_bar(aes(health, contact, fill = Condition),
           stat = "identity",
           color = "black",
           position = position_dodge(.9),
           alpha = .8) +
  geom_errorbar(aes(ymin = contact - ci, ymax = contact + ci),
                width = .25,
                position = position_dodge(.9),
                alpha = .6) +
  theme_classic(base_size = 20) + 
  scale_fill_manual(values = c("gray70", "turquoise4")) + 
  scale_color_manual(values = c("gray70", "turquoise4")) +
  facet_wrap(~ Source, labeller = labeller(Source = sourceLabs)) + # this is that simple trick
  xlab("") +
  ylab("Willingness to Have Contact with Targets") +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels = c("Diseased Targets", "Healthy Targets")) 

#  ggsave("contact means plot (with points).jpg", device = "jpeg", units = "cm")


# #  run additional t-tests against midpoint (0) for aversion to any face for crowded and control conditions
# crowded/healthy 
crowdHealthy <- filter(cdLongA, Condition == "Crowd" & health == "Healthy")

# test
t.test(crowdHealthy$contact, mu = 0)
cohen.d(crowdHealthy$contact, f = rep(0, each = nrow(crowdHealthy)))


# crowded/diseased
crowdDisease <- filter(cdLongA, Condition == "Crowd" & health == "Diseased")

# test
t.test(crowdDisease$contact, mu = 0)
cohen.d(crowdDisease$contact, f = rep(0, each = nrow(crowdDisease)))


# control/healthy
contrHealthy <- filter(cdLongA, Condition == "Control" & health == "Healthy")

# test
t.test(contrHealthy$contact, mu = 0)
cohen.d(contrHealthy$contact, f = rep(0, each = nrow(contrHealthy)))


# control/diseased
contrDisease <- filter(cdLongA, Condition == "Control" & health == "Diseased")

# test
t.test(contrDisease$contact, mu = 0)
cohen.d(contrDisease$contact, f = rep(0, each = nrow(contrDisease)))


# visualize these
zeroSum <- Rmisc::summarySE(cdLongA,
                            measurevar = "contact",
                            groupvars = c("health", "Condition"),
                            na.rm = T)
zeroSum

ggplot(zeroSum, aes(health, contact, fill = Condition)) +
  geom_point(data = cdLongA,
             aes(health, contact, color = Condition),
             position = position_jitterdodge(.15, .15, .9),
             alpha = .5, size = .75) +
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(.9),
           alpha = .8) +
  geom_errorbar(aes(ymin = contact - ci, ymax = contact + ci),
                width = .25,
                position = position_dodge(.9),
                alpha = .6) +
  theme_classic(base_size = 20) +
  scale_fill_manual(values = c("gray70", "turquoise4")) +
  scale_color_manual(values = c("gray70", "turquoise4")) +
  xlab("") +
  ylab("Willingness to Have Contact with Targets") +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels = c("Diseased Targets", "Healthy Targets")) 

#  ggsave("zero comp plot.jpg", device = "jpeg", units = "cm")
