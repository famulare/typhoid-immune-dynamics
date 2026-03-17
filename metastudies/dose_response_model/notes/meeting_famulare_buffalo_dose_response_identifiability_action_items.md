# Action Items: Dose-Response Identifiability Discussion

**Source**: `meeting_famulare_buffalo_dose_response_identifiability.md`

---

## Vince

1. Relax the currently fixed field-dose parameter in the existing model and verify that the expected exposure-rate/field-dose non-identifiability appears in the posterior.
2. If the model appears more identifiable than expected, investigate what information, parameterization, or geometry is creating that behavior.
3. Extend the data-generation and identifiability pipeline to include fever as a separate outcome channel rather than treating the model as infection-only.
4. Use the existing ABM structure, where possible, to expose the infection-versus-fever pattern already present in the system.
5. Build the Stan-side identifiability stack far enough to demonstrate where added outcome channels and added biological parameters create new non-identifiabilities.
6. Draft a written plan based on the transcript and send it back to Mike for review.
7. Share initial results by the end of the day or early the next morning, ahead of the 11:30 meeting.

---

## Mike

1. Continue reasoning carefully through the available data to determine how, and to what extent, `N50`, `alpha_DR`, `gamma`, and the infection-versus-fever split can actually be identified.
2. Clarify which assumptions are scientifically defensible for constraining the infection and fever channels, especially whether `N50` should be shared and which parameters are allowed to differ.
3. Trace which elements of the current assumptions are inherited from prior polio and vaccinology work and assess how portable those assumptions are to typhoid.
4. Map the real datasets that could provide the orthogonal information needed beyond incidence alone.
5. Distinguish clearly between what is biologically identifiable in an ideal clean-data setting and what is realistically identifiable with the available typhoid data.
6. Prepare guidance for the next discussion on how far the model can be pushed before identifiability depends mainly on priors or algebraic constraints.

---

## Shared Understanding To Carry Forward

- Incidence-style field data alone identify exposure rate naturally, but do not separately identify field dose once that parameter is introduced.
- Breaking that structural non-identifiability requires orthogonal information, not just more incidence data.
- The most promising additional information streams are multi-dose CHIM data and outcome splits between infection and fever across immune states.
- Any expanded model will likely require explicit scientific constraints, not just more flexible computation.
