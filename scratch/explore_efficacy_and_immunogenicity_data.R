# efficacy and immunogenicity exploration scratchfile

library(tidyverse)

# efficacy data
d_eff = readxl::read_excel('data/typhoid_vaccine_study_data.xlsx',sheet = 'efficacy')  |>
  mutate(estimate=as.numeric(estimate),
         lower=as.numeric(lower),
         upper=as.numeric(upper),
         age_min_years = as.numeric(age_min_years),
         age_max_years = as.numeric(age_max_years)) |>
  mutate(timepoint_months_authors_mean=1/2*(timepoint_months_authors_min + timepoint_months_authors_max)) |>
  mutate(age_mean = 1/2*(age_min_years + age_max_years)) |>
  mutate(age_label = reorder(factor(paste(age_min_years,'-',age_max_years,sep='')),age_mean)) |>
  mutate(age_label_coarse=cut(age_mean,breaks=c(0,2.5,6.35,10,15,40),
                              labels=c('toddler','younger child','children','early teen','late teen'),ordered_result=TRUE))
  
  
names(d_eff)

# plot total individual VE for bacteriologically-confirmed typhoid fever vs time
plot_efficacy_df = d_eff |> 
  filter(endpoint == 'total_individual_VE') |>
  filter(outcome %in% c('blood_culture_confirmed_typhoid_fever')) |>
  filter(vaccine == 'Typbar-TCV') |>
  arrange(study,age_label,timepoint_months_authors_mean)

names(plot_efficacy_df)

plot_dat = plot_efficacy_df |> mutate(group=interaction(age_label,design_group)) |>
  select(dotname,age_label_coarse,group,timepoint_months_authors_min,
         timepoint_months_authors_max,estimate,estimator,endpoint_duration,age_cohort,design_group,
         lower,upper,timepoint_months_authors_mean)

plot_dat = plot_dat |>
  rbind(plot_dat |> group_by(group) |>
          mutate(timepoint_months_authors_min = timepoint_months_authors_max,
                 timepoint_months_authors_mean=NaN))

ggplot(plot_dat, aes(shape=dotname,color=age_label_coarse,fill=age_label_coarse,group=group)) +
  geom_ribbon(aes(x=timepoint_months_authors_min/12,ymin=lower,ymax=upper),color=NA,alpha=0.1) +
  geom_step(aes(x=timepoint_months_authors_min/12,y=estimate,linetype=estimator)) +
  geom_point(aes(x=timepoint_months_authors_mean/12,y=estimate)) +
  facet_grid('endpoint_duration ~ age_cohort') +
  theme_bw()  +
  scale_color_brewer(palette='Paired') +
  scale_fill_brewer(palette='Paired') +
  xlab('years since vaccination') +
  ylab('VE estimate')
ggsave('scratch/figures/efficacy_time_age.png',units='in',width=7, height=5)



ggplot(plot_dat, aes(color=dotname,fill=dotname,group=group,shape = age_label_coarse)) +
  geom_ribbon(aes(x=timepoint_months_authors_min/12,ymin=lower,ymax=upper),color=NA,alpha=0.1) +
  geom_step(aes(x=timepoint_months_authors_min/12,y=estimate,linetype=estimator)) +
  geom_point(aes(x=timepoint_months_authors_mean/12,y=estimate)) +
  facet_grid('endpoint_duration ~ age_cohort') +
  theme_bw() +
  xlab('years since vaccination') +
  ylab('VE estimate')
ggsave('scratch/figures/efficacy_time_location.png',units='in',width=7, height=5)


# immunogenicitiy data
d_imm = readxl::read_excel('data/typhoid_vaccine_study_data.xlsx',sheet = 'immunogenicity') |>
  pivot_longer(
    cols = matches("^(treatment|control)_day"),    # select columns starting with "treatment" or "control" and "day"
    names_to = c("cohort", "day", "measure"),       # create three new columns from parts of the name
    names_pattern = "(treatment|control)_day(\\d+)_(.*)",  # regex to extract the cohort, day, and measure
    values_to = "value",                             # the data goes into a new column named 'value'
    values_transform = list(value = as.character)
  ) |>
  mutate(day = as.numeric(day)) |>
  pivot_wider(
    names_from = measure,
    values_from = value
  ) |>
  drop_na(N) |>
  mutate(across(c("age_min_years","age_max_years","N", "seropositive_percent", "fold_rise_median","fold_rise_lower",
                  "fold_rise_upper","seroconversion_percent"), as.numeric)) |>
  mutate(elisa_U_per_ml_median_numeric = if_else(elisa_U_per_ml_median=='<7.5',3.75,if_else(elisa_U_per_ml_median=='<3.0',1.5,as.numeric(elisa_U_per_ml_median)))) |>
  mutate(elisa_U_per_ml_lower_numeric = if_else(elisa_U_per_ml_lower=='<7.5',3.75,if_else(elisa_U_per_ml_lower=='<3.0',1.5,as.numeric(elisa_U_per_ml_lower)))) |>
  mutate(elisa_U_per_ml_upper_numeric = if_else(elisa_U_per_ml_upper=='<7.5',3.75,if_else(elisa_U_per_ml_upper=='<3.0',1.5,as.numeric(elisa_U_per_ml_upper)))) |>
  mutate(age_mean = 1/2*(age_min_years + age_max_years)) |>
  mutate(age_label = reorder(factor(paste(age_min_years,'-',age_max_years,sep='')),age_mean))

names(d_imm)

# plot total individual VE for bacteriologically-confirmed typhoid fever vs time
plot_dat = d_imm |>
  filter(cohort == 'treatment') 
 
ggplot(plot_dat, aes(shape=dotname,color=age_label,fill=age_label,group=age_label)) +
  geom_ribbon(aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric),color=NA,alpha=0.1) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  facet_grid('antigen ~ vaccine') +
  theme_bw() +
  scale_y_continuous(trans='log10')

ggplot(plot_dat, aes(shape=dotname ,color=dotname,fill=dotname,group=age_label)) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  facet_grid('antigen ~ vaccine') +
  theme_bw() +
  scale_y_continuous(trans='log10')


