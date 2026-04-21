# Scientific Narrative: Dose-Response Identifiability Discussion

**Source**: `meeting_famulare_buffalo_dose_response_identifiability.md`  
**Participants**: Mike Famulare, Vince Buffalo  
**Purpose**: Reconstruct the scientific logic of the conversation in a form that can be read without having attended the meeting.

---

## Core Question

The discussion centered on a basic identifiability problem in the typhoid dose-response framework: once the model includes both an ecological exposure rate and an ecological field dose, what information is available to distinguish them, and what additional data or assumptions are needed to identify the biology and the ecology separately?

---

## 1. What the current incidence-style data can identify

The starting point was the WHO-style scenario data, treated as if each scenario were a real setting with a characteristic average incidence level. In a standard epidemiologic model without an explicit dose-response layer, these data naturally support one ecological degree of freedom per setting: an exposure rate. With three settings, for example, one would fit three exposure-rate parameters while keeping biological parameters such as immune kinetics fixed across settings.

That setup is straightforward and is essentially equivalent in spirit to a more conventional transmission model in which settings differ only in transmission intensity.

---

## 2. Why adding field dose creates a strong non-identifiability

The dose-response framework introduces a second ecological parameter: field dose. If the only information coming from the field data is incidence plus antibody dynamics, then exposure rate and field dose become strongly non-identifiable. For any setting, there is a continuum of exposure-rate/field-dose combinations that produce the same effective infection probability and therefore the same likelihood.

This is not a weak identifiability issue that might be resolved by better sampling or tighter computation. It is a structural problem: the likelihood contains a ridge. In practice, Stan should discover that ridge as a strong posterior dependence between exposure rate and field dose.

The immediate implication is that field dose cannot simply be freed and estimated from these data alone. Additional information is required that is not just more of the same incidence signal.

---

## 3. The kind of information needed to break the ridge

The key scientific point was that the model does not merely need more data; it needs orthogonal information. Specifically, it needs information that distinguishes biological dose-response behavior from ecological exposure frequency.

Two sources of such information were discussed.

### 3.1 Multi-dose CHIM data as the cleanest source of biological identification

The cleanest possible source would be a modern multi-dose controlled human infection model in clearly typhoid-naive adults, spanning a wide range of doses and ideally delivered in a realistic food vehicle rather than only a bicarbonate preparation. In that idealized setting, the challenge data could identify the naive-host dose-response curve directly, especially parameters such as `N50` and `alpha_DR`.

The reason is that multiple challenge doses over a wide range reveal the shape of the dose-response curve itself. Under the working assumption that basic host biology is portable across settings while ecological conditions are not, those biologically interpretable parameters could be estimated independently of transmission ecology.

This ideal dataset does not exist in the desired form, but it serves as the conceptual benchmark for what "orthogonal information" would look like.

### 3.2 Infection-versus-fever data as a second orthogonal signal

Even if naive-host dose-response parameters were identified from CHIM studies, that would still not separate field dose from exposure rate in real settings. A second layer of information is needed.

The promising source discussed here is the distinction between infection and fever. Both CHIM studies and surveillance suggest that some individuals show evidence of infection without fever, and that the infection-to-fever relationship changes with immunity. That divergence creates leverage because it provides information on two related but distinct outcomes, both modulated by immune status.

In the meeting, this outcome split was treated as potentially informative about:

- the immune-scaling parameter `gamma`, because immunity changes the gap between infection and fever
- additional dose-response structure, because infection and fever may not share exactly the same dose-response curve

The broad idea is that once the model admits both infection and fever as outcomes, and once their divergence across immunity levels is informed by challenge data, the field data may contain enough structure to estimate both ecological parameters per setting.

---

## 4. The modeling assumptions under discussion

A major part of the conversation concerned which algebraic or biological constraints might be needed to make the expanded model identifiable.

One recurring assumption was that `N50` may be shared across infection and fever, while other parameters differ between the two outcome channels. Mike noted that in prior work, especially ideas imported from polio and vaccinology, he had assumed a common `N50` but allowed `alpha_DR` to differ between infection and fever. He also described `gamma` as a key immune-dependent quantity informed by the infection-fever divergence, and later indicated that infection and fever might need separate `gamma` parameters as well.

The important takeaway is not that these assumptions were fully settled. They were not. The meeting instead established that:

- an unconstrained infection/fever extension adds multiple degrees of freedom quickly
- without rich challenge data, several of those degrees of freedom will remain non-identifiable
- identifiability may therefore depend on explicit structural assumptions, informative priors, or algebraic constraints, not just more computation

This makes the model-building problem partly scientific and partly statistical: the team needs to decide which constraints are biologically defensible before asking Stan to estimate the full system.

---

## 5. The intended information flow for the full framework

The meeting outlined an ideal sequence for inference:

1. Use clean challenge data to identify portable biological parameters such as the naive dose-response curve and immunity-dependent modifiers.
2. Use infection-versus-fever information, especially across immune states, to identify how immunity shifts outcomes and how the two outcome channels relate.
3. Only after those biological pieces are reasonably pinned down, fit ecological parameters such as field dose and exposure rate in each real-world setting.

Under this logic, ecology should be learned only after biology is sufficiently constrained. Otherwise the model tries to estimate too many interacting quantities from incidence-like data alone and collapses into structural non-identifiability.

---

## 6. Why typhoid is harder than the original polio use case

Mike emphasized that this framework is conceptually cleaner in polio, where oral vaccine challenge data act like well-characterized challenge studies in well-defined cohorts. That level of clarity creates a more effective "flywheel" between challenge data, biological parameterization, and field inference.

Typhoid is harder because the corresponding data are messier. Cohorts are less cleanly characterized, outcomes are less standardized, and the available challenge evidence is not as rich as one would want for direct identification of the biological parameters. This is why the conceptual model can be stated clearly while the rigorous Stan implementation remains difficult.

---

## 7. The immediate computational implication

Vince proposed an important near-term test: simply relax the currently fixed field-dose parameter in the present model and check whether the expected non-identifiability appears. The expectation is that it should. If it does not, that would be scientifically interesting because it would imply that some other feature of the model or likelihood is supplying unexpected information or imposing an accidental constraint.

This is an important diagnostic step because it separates a conceptual identifiability claim from actual posterior behavior in the implemented model.

---

## 8. Where the discussion landed

By the end of the meeting, the group had a clearer division between conceptual and practical next steps.

Conceptually, the meeting clarified that the hard problem is not merely fitting a richer model; it is specifying which data streams and which assumptions are capable of separating:

- portable biological dose-response parameters
- immune modulation of outcomes
- setting-specific ecological exposure rate
- setting-specific field dose

Practically, the next stage is to make the expected non-identifiabilities visible in the computational pipeline and to map exactly which real data sources could justify the constraints needed to recover identifiability.
