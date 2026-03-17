# Meeting: Dose Response Model Identifiability Discussion

**Participants:** Mike Famulare, Vince Buffalo

---

**Vince Buffalo:** So what does the data look like in those cases? That's just more of the — it's not a challenge study, it's like which of the modules would that type of data fit into?

**Mike Famulare:** Well, let's — yeah, there's two different pieces. So I'm gonna start with what's closest to what you've been doing, and I think what I need to do is like, in order to get the ideas out coherently, I'm kind of gonna do my LLM autoregressive generation of the text right now.

So you know, the first type of data that's informative is very much like the WHO scenario setup. So you know, let's continue to ignore that those are not real places but are composites of places. Let's treat them as real places, right? For the sake of discussion. So you know, those data come with an incidence level, right? There's, you know, even ignoring age for a moment, each of those scenarios has an average incidence level.

And that is exactly one of your inputs, right? And so from those scenarios, you don't have anything to power the negative binomial dispersion. But that's OK, whatever — treat it as Poisson. That's a nuisance parameter. Anyway, you have an overall level of incidence, right.

So the two things — like right now, if you were not operating in this dose response framework, if you were just operating in your sort of standard epi model framework like Kira and them are doing right now, you would have just a different exposure rate for each setting, right? That's the natural degree of freedom. It is the natural ecological single parameter if you only have a model that admits one ecological parameter. The parameter for transmission — it's the exposure rate, right?

And so you trivially at this point, and I think this is equivalent, you know, intellectually equivalent to what Kira has already done — you would have 3 settings, you'd fit 3 exposure rates, and everything else that you think is not setting-specific like T_to_K and kappa and so on, you would keep fixed across settings.

And so that's what they've kind of already done. Extra bells and whistles because the model has extra features, but that's sort of where they come from, right? Does this make sense?

**Vince Buffalo:** Mm-hmm.

**Mike Famulare:** Cool. So the dose response framework says — without bringing in any other data — the dose response framework introduces a second ecological parameter, which is the field dose. And with no other information besides incidence and antibody dynamics, now you've got no ecological identifiability, right? Because every given site has a hyperbola between exposure rate and field dose that yields the same probabilities and so yields exactly the same results.

And so if you tried to just free up field dose, now with three sites or whatever, each site would have — you would co-vary field dose and — you have 6 parameters for three sites, right? So it's obviously unidentifiable, and your Stan will find a, you know, whatever — some non-transformed hyperbola-ish or transform line thingy, right — a ridge in likelihood space.

So like you did with your synthetic challenge data, one has to find data that can independently distinguish those — or I should be careful here. Let me be careful. One needs *information* — I shouldn't say data — that can distinguish field dose from exposure rate.

**Vince Buffalo:** I mean these — so those, the CFUs, exposure rate is like a hazard, like a probability or a rate. So — and they're really, I mean, these are capturing the same thing as you're saying. This isn't like a weak identifiability problem. It's like this thing — these are scalar or functional differences.

**Mike Famulare:** This is a strong identifiability problem. Yeah, they are — strict. Without extra information at this point, they are — there is a, like you already figured out, there's an arithmetic relationship between them that yields the same likelihood. Right? So there's like there is a function that has infinity — infinite pairs equal the same likelihood — that is just a function-inversion-of relationship between exposure rate and field dose. So this is a strong non-identifiability problem at this point.

All right, so you have to get some information in that breaks that identifiability problem. It has to bring in information that's not incidence, at least not incidence in the same way — setting-specific, just like how many infections a month.

So if there were no complications to the data — which, so I'm tabling complications for a moment, which are real, and this is why Kira had the questions that she had for me — you know, if for example it was the challenge, so the CHIM studies, particularly CHIM studies that provide multiple dose challenge, that is the most obvious source of orthogonal information.

So like, under the factoring we talked about before — that we're going to assume individual biological things are portable, and ecological things aren't — you would have, let's say you imagine that the Hornick study was, you know — I mean, to be honest, imagine that the Hornick study was an Oxford study in 2016, but they did it with the thoroughness that Hornick did. That would be the ideal. I would love to have that study. That doesn't exist. I would love to have definitely-not-exposed adults, you know, that definitely got 9 different doses or seven different doses across 8 orders of magnitude, right? And delivered it in a food vehicle like milk, not some bicarbonate solution.

Like ideally I would have milk typhoid from 10^1 to 10^9, and I would have it in adults from England that have definitely never had typhoid before, right? Like if I had that, that would be great orthogonal information to completely determine N50 and alpha_DR. Because that curve over that range is exactly the kind of sigmoid that's not actually an exponential 1-e^(-x) — it's exactly the beta-Poisson functional form, right? So and by assumption, adults with no exposure history are a portable population of adults. You wouldn't even need a joint likelihood for this. If it was truly that simple, you would just fit out N50 and alpha_DR, and then you would have those completely identified separately, right? So that's the dream.

Now then to come back to the incidence thing, you still don't have orthogonal information about exposure and dose.

The second piece of orthogonal information that I was using in the prototyping is a distinction between fever and infection. Because a known feature from the CHIM studies in typhoid and clinical surveillance, although less so, is that you can have evidence for infection — whether it be seroconversion events or stool positivity or blood culture positivity — and not fever, and that as a function of immunity.

Those ratios of infection to fever do change, right? So in your data space, you have something that's informative about two different outcomes and depends on immunity. So that gives you — with two different outcomes and a dependence on immunity, you got two more pieces of information. There are two more degrees of freedom that you can pin down.

One of those degrees of freedom should be highly informative about gamma, because we know it is an immune-dependent phenomenon. The divergence between infection and fever is immune-dependent, and that'll be directly informative about gamma.

The other degree of freedom is of course that some combination of N50 and alpha_DR should be indexed by infection versus fever, right? Because they are different manifestations of the same condition. I made an assumption that we can talk through, although I don't want to do it right now because I will lose the rest of these threads.

We will assume that alpha_DR is the thing that — of the constraint surface between N50 and alpha_DR for fever versus infection — I assume that N50 was the same for both and that alpha_DR was different. And so that's the extra degree. So I've added another parameter because I took alpha_DR from one to two, and I've got gamma, and I've got now the dose response for — we've got dose response for naive people, dose response for fever, you've got fever versus infection as a function of immunity. So that brings in non-immune people, and you've got the range of both dose and a range over immunity.

And so you can kind of — I hopefully can kind of see how like, OK, I've got 4 parameters, I've got two variable-rate data pieces. Each piece of data is roughly informative about two of those 4 parameters. Or in some projection of four parameters, each of which is informed by two pieces of data. I can nail down under this assumption that I have this good of data — I can nail down these assumptions about N50, alpha, and gamma before I ever start fitting a transmission ecology.

And then I still have field dose and exposure rate, so I still have two parameters for every ecology, but I've allowed into my state space fever and infection. And now for every ecology I have two measurements. And so now I can fit both — I can fit those separate from exposure, and the information about that comes from having two measurements and their divergence as a function of immunity, which I've already pinned down by the CHIM studies.

So like, this is the ideal information flow. And in polio, where this model was originally developed, we have the cleanliness of these cohorts being this well defined because we can use OPV — we use the oral vaccine — and so we have, you know, vaccine-quality CHIM studies. They're very well characterized background cohort populations.

And so that's the flywheel. And then the hard part, and the reason Kira has valid doubts, and the reason I have valid doubts — the difference between how I solve these things hand-fit and how we would really solve them if we're really rigorous and take the seriousness of Stan — is how do we deal with the things that aren't that clean.

Yeah, so I just threw a lot at you. I'm glad we have the transcript. How you doing?

**Vince Buffalo:** Yeah, I mean, I think it would be helpful to translate this into current gaps in the sort of data generation pipeline. I mean I think, you know, given that it's been fixing the field dose response right to fit the stuff initially, that model — relaxing that, seeing it behave poorly — those feel like next steps to me. But thinking about it in the way that you just mentioned, like having something that generates ecological-like data, something that separates infection from fever, having that in the data generating process feels to me like the next step.

I think we can take the transcript and I can give it to my Claude and then sort of say like, this is what we're thinking, and then I can send back a plan or something to you and you can sort of see how that looks. I think before I do that though, it is useful to relax that one parameter. Sometimes, or just make it unpinned and just make sure that it's unidentifiable. You never know — sometimes the posterior geometry's weird and it just latches on to other information. But I agree there shouldn't conceptually be any information that isn't there, but you never know.

Like, yeah, in the Klebsiella case, I will just say — there it was interesting because the data was fitting very poorly, like we couldn't match the observed data in a few cases. And that's a case where it's lots of different studies, different confounders, different sort of conditioning on population, et cetera. It's a bunch of sources of ascertainment bias, blah blah blah. But Dan had added in one extra parameter, and so the Stan did this fit and it worked pretty well. Like, it was much better than any of the ABM-based calibration approaches.

Stan — or Dan — added in one extra parameter and my Stan (or my Claude) caught that, and it was just like, all right, let's add this parameter in and it refit the model. The number of parameters went from 7 to 8 when we added in that new parameter, but then the number of effective parameters went from 7 to 3, right?

**Mike Famulare:** I remember hearing you say this, yeah.

**Vince Buffalo:** Yeah. And it was just — it was a truly incredible moment of model building in the Bayesian world, because it's just like, oh, now we actually don't need to strain all 7 or 8 dimensions of the parameter space to get it to fit the data. Really, there's just some simpler manifold that it's able to move across.

So that was, I think — I haven't done a lot of model comparison. I've done the across-model comparison with different data, but that's not really a model comparison activity because you can't compare things with different data. But now that we refine the model, I think hopefully by the end of today, I can push on that and share some stuff. What time is the meeting tomorrow?

**Mike Famulare:** 11:30 or something like that.

**Vince Buffalo:** Yeah, 11:30. Yeah, I should be able to have some results by the end of today and early tomorrow I can share. But does that sound like a good plan?

**Mike Famulare:** Yeah, let me just repeat a little bit of it back to you so I make sure I understand it, cause that tells me the pieces I need to take on for tomorrow as well.

OK, so your plan as I understood it is something like — the first thing to do is just immediately relax the field dose parameter. And expect to find the unidentifiability naturally. But if not, that would be interesting to make sure we figure out why.

Then the second thing is in the data generating process, add fever.

**Vince Buffalo:** Yeah, I figured it was already there.

**Mike Famulare:** Yeah — you go back to the ABM that you have, that should already be there. So the pattern for that should already be there.

And then — so when I say "add it" in the context of you using Stan in this model-building approach that you're using, which I think is excellent and different than the way I tend to do it — because I do it in my head more, and I like the computational rigor you're bringing to this a lot — in principle, adding a fever dose response is going to add a new N50, a new alpha_DR, and a new gamma, right? So we're adding one channel with three parameters.

And so in your identifiability pipeline, there's a mix of computational identifiability constraints — like how would I constrain these things to re-identify it? — and also the scientifically — given that in principle there's a constraint surface here, it's a three-parameter thing, it's completely under-constrained mathematically. The equivalent — in Stan it could be literally just priors that are the constraints, or go one step further than just priors to algebraic assumptions.

**Vince Buffalo:** *(had to step away — dog needed to go outside)*

**Mike Famulare:** *(continuing as monologue for the transcript)*

And yeah, so then we're talking about multiple degrees of freedom. So I don't fully remember my reasoning for why I picked what I picked, but I made the assumption — I think I did this back — actually, I know I did this. This is laundered from polio, where the data support the decision I made. And it's laundered from vaccinology. I do actually know how to prove this.

I make an assumption that N50 is the same across infection and fever, but alpha for infection and fever are allowed to be different. And then I also make an assumption that gamma for infection and fever are allowed to be different. And without rich CHIM data, those will be multiple non-identifiabilities.

And so I think what I hear for your part — for today and tomorrow — is to sort of build up that Stan stack that gets us to that point where those non-identifiabilities are clearly demonstrated.

And my part of this for the afternoon is to continue to reason through the data carefully, to really plan out like how and to what extent can we identify them.
