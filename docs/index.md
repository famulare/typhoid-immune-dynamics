# Project overview: typhoid immune dynamics

## What

We are building a new typhoid immunity model with a continuous, time-varying correlate of protection and dose-dependent susceptibility.

## Why

Because typhoid immunity is interesting! The new Typhoid Conjugate Vaccine (TCV) seems to perform differently in different settings. Everywhere it's been trialed, the observed efficacy against blood-culture confirmed typhoid is around \~80% for the first two years, but while there is little evidence of waning in some settings — not yet reported in [Nepal](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(21)00346-6/fulltext) and statistically insignificant in [Malawi](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850983/) — there is clear waning in [Bangladesh](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(24)01494-6/fulltext) and both lower observed efficacy and faster and city-dependent waning in [Pakistan](https://www.medrxiv.org/content/10.1101/2024.08.30.24312839v1.full). In addition, there is some observed variation of waning by age, with younger children waning faster.

These strong differences in vaccine effectiveness by setting and age indicate the need for bespoke vaccination strategies. To rationally design those strategies across settings — from hyper-endemic to solely at risk from importation — we need a model capable of explaining the observations and extrapolating from them.

## How

[Prior work by IDM on poliovirus](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2002468) — another predominantly enteric pathogen transmitted via the fecal-oral route — found that epidemiological measures of immunity depend on both individual-level immune response and the typical size of pathogen doses acquired via transmission, and we built a model to capture that dynamic. Following the demonstration that a [similar model framework was useful for COVID](https://pubmed.ncbi.nlm.nih.gov/36976678/) as well, we've developed a general formulation of the model and will be fitting it to typhoid for the first time in this project.

We'll explain the model as we develop it using this blog and documentation. Until then, here are preview [slides describing the approach](assets/Better-defaults-in-acquired-immunity-modeling.pdf). You can see that it is in the flavor of a PK/PD model common for drug development, but with some features I haven't seen implemented by others. For a technical introduction to PK/PD modeling, see these [notes by Henrik Madsen](https://www.henrikmadsen.org/wp-content/uploads/2014/10/Report_Peer_reviewed_-_2008_-_Introduction_to_PK_PD_modelling_-_with_focus_on_PK_and_Stochastic_differential_equations.pdf)).

## When

We hope to have useful results well before the WHO SAGE meeting on TCV schedule recommendations in October 2025. Initial deadline to de-risk the whole project is April 2025.

## Where

Developed at the [Institute for Disease Modeling (IDM)](https://www.idmod.org/), a research institute within the [Gates Foundation](https://www.gatesfoundation.org/our-work#program_strategies)'s Global Health Division.

The Github repo is [famulare/typhoid-immune-dynamics](https://github.com/famulare/typhoid-immune-dynamics).

## Who

Technical project lead: [Mike Famulare](https://scholar.google.com/citations?user=TPWwr18AAAAJ&hl=en).

Technical collaboration with: [Kyra Grantz](https://scholar.google.com/citations?user=pDS-Fk8AAAAJ&hl=en) and [Alicia Kraay](https://scholar.google.com/citations?user=Qc2kca0AAAAJ&hl=en&oi=ao).

IDM Typhoid research and strategy lead: [Jillian Gauld](https://pubmed.ncbi.nlm.nih.gov/?term=%28%28Gauld%2C+Jillian%5BAuthor%5D%29+OR+%28Gauld%2C+JS%5BAuthor%5D%29%29+NOT+%28Gauld%2C+JW%5BAuthor%5D%29&sort=pubdate).
