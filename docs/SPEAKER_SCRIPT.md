# RO-2026-008 | Speaker Script
## BU Health Data Science & AI Showcase, 15 September 2026, Hiebert Lounge

Four versions, layered. Learn the 3-minute one properly; the others are
compression and expansion of the same spine.

**Spine, memorise this and everything else follows:**

> Decontamination asks whether a microbe is really there. Nobody was asking
> whether the *method* could tell if it mattered. So I built data with nothing
> in it. The standard tests passed anyway. Then I checked a bacterium where we
> already know the answer, and the method got the direction backwards.

**Delivery notes**
- Stand to the *left* of the poster so your body does not block column 2, panel 3. That is the panel that sells it.
- Do not read the poster aloud. They can read. Talk about what is not written.
- Point with an open hand, not a finger.
- When you hit the H. pylori panel, slow down. That is the moment.
- If someone stops you mid-sentence with a question, take it. A question means they are engaged; finishing your paragraph does not.

---

# VERSION 1 — Thirty seconds (walk-by)

Use when someone slows but does not commit. One breath, then let them leave or stay.

> "Short version: the tumour microbiome field found bacteria inside cancers and
> started asking which ones cause it. I asked a step back, whether our methods
> could even tell. So I built a fake dataset with no biology in it at all. The
> two tests this field runs passed on it, ten times out of ten.
>
> Then I checked *H. pylori*, which we know causes stomach cancer. In real
> data it's eight times *less* abundant in the tumour. For the one organism
> where we know the truth, the standard approach gives you the opposite
> answer."

Then stop. If they lean in, go to Version 2.

---

# VERSION 2 — Three minutes (standard walkthrough)

## Opening, 20 seconds

> "Thanks for stopping. Two-minute version, stop me anywhere.
>
> In 2024 a major pan-cancer tumour microbiome paper was retracted, because
> what looked like biological signal turned out to be human sequencing reads
> misclassified as bacteria, plus batch effects. The field responded well: it
> built decontamination pipelines and got much better at asking *is this
> microbe really here*.
>
> But that leaves a second question untouched. If it *is* here, is it causing
> anything, and which way does the arrow point?"

## Panel 1, the gap, 20 seconds

*[gesture, column 1 top]*

> "Almost all of this work compares tumour tissue to healthy tissue next door
> and asks which bacteria differ. That's a photograph. Causation is a movie.
> You can't tell from one frame whether the bacterium arrived and caused
> trouble, or arrived because there was already trouble."

## Panel 2, the audit, 40 seconds

*[move to Figure 1]*

> "So before touching real patients, I asked whether we'd even notice being
> wrong.
>
> I built two synthetic datasets. Same size, same shape, same code. One has a
> real pattern I planted. The other has nothing, by construction. And in the
> empty one I deliberately tangled which lab processed a sample with whether it
> was labelled cancer, because that's what happens in the real world.
>
> Then I ran the two tests this literature reports most: comparison against a
> no-information rate, and label permutation. Both passed on the dataset
> containing nothing. Ten seeds out of ten.
>
> Only two checks caught it, and neither is commonly reported: comparing
> against a model built from the confounder alone, and training on some labs
> and testing on a lab the model has never seen."

*[tap the 2.5x annotation]*

> "Batch confounding on its own drove accuracy to two and a half times chance,
> from data with no biology in it."

## Panel 3, the inversion, 45 seconds — SLOW DOWN

*[move to column 2, Figure 4b]*

> "Then I wanted a real-world test, not a simulation. There's one organism
> where we already know the answer. *Helicobacter pylori* causes stomach
> cancer. That's not controversial, it's a recognised carcinogen and treating
> it prevents the disease.
>
> So: in a careful, decontaminated dataset, is it enriched in the tumour?
>
> It's about eight-fold *depleted*. Thirty-nine patient-matched pairs, same
> patient, tumour versus their own adjacent normal tissue, p equals four times
> ten to the minus four."

*[pause two beats]*

> "And that makes complete sense once you know the biology. *H. pylori* needs
> a highly acidic stomach lining to survive. The cancer process destroys
> exactly that lining, through atrophy and intestinal metaplasia. So the
> bacterium starts the fire, and the fire burns down the house it was living
> in.
>
> It's the arsonist who isn't at the scene when the fire trucks arrive. If
> your method is 'look who's standing near the fire', you arrest the
> firefighters."

## Panel 5 and 6, the fix, 45 seconds

*[move to column 3]*

> "So what do you do instead? You need something the disease cannot have
> caused. Your genes work. You had them at conception, decades before any
> tumour, so a tumour can't change them. And some people carry variants that
> nudge them toward hosting more of one bacterium. That's a coin flip nature
> performed at birth.
>
> I did that for 211 gut bacteria against colorectal cancer, then repeated it
> in a completely separate cohort of a hundred thousand cases.
>
> Honest answer: no single bacterium held up. Every promising lead from the
> first cohort vanished in the second. I'm not going to tell you I found the
> colon cancer bug.
>
> But something did survive. Comparing the two independent cohorts, the
> bacteria point the *same direction* about 61 percent of the time, where
> chance would give you 50. So the signal is real. It's just spread thinly
> across many organisms, none individually strong enough to find at current
> sample sizes. A choir, not a soloist."

## Close, 15 seconds

> "Two practical takeaways. If you see a headline about the bacterium behind
> some cancer, and the design was tumour versus healthy tissue, that finding
> may be real but the design cannot tell you which way the arrow points. And
> when a study reports nothing, ask how big an effect they could have detected.
> Mine was 0.8 log units. Saying you found nothing without that is like saying
> you searched the house without mentioning you never opened any doors.
>
> What did I lose you on?"

---

# VERSION 3 — Ten minutes (lecture)

Use if someone sits down with you, or for a seminar slot. Same spine, with the
mechanism spelled out. Plain language throughout; assume no microbiology and no
statistics beyond an average.

## 1. Why this problem exists at all, 90 seconds

> "Let me start with why anyone thought there were bacteria in tumours.
>
> When you sequence a tumour, you get hundreds of millions of short DNA
> fragments. Most are human. Some are not. If you take the leftovers and match
> them against a database of bacterial genomes, you get hits. In 2020 people
> did this across many cancers and found what looked like a distinctive
> bacterial signature for each cancer type. It was exciting: a blood test for
> cancer based on microbes.
>
> The problem is that 'leftover fragment that matched a bacterial genome' is
> not the same as 'bacterium was present'. Human DNA can match bacterial
> databases by chance, especially in repetitive regions. Laboratory reagents
> carry bacterial DNA. And crucially, samples from different cancers were
> often processed at different centres, so anything that differs between
> centres looks like it differs between cancers.
>
> In 2024 the main paper was retracted. The field's response was correct: build
> better decontamination. Filter the reagent contaminants, correct for batch,
> use stricter matching. Good work, and I use its output.
>
> But decontamination answers one question. It tells you the microbe is really
> there. It says nothing about whether it does anything."

## 2. The two questions, and why the second is harder, 90 seconds

> "Here's the distinction, and it's the whole talk.
>
> Question one: is this bacterium present? That's a measurement problem. You
> solve it with better filtering and better controls. The field has largely
> solved it.
>
> Question two: does it matter, and which direction? That's a *causal* problem,
> and no amount of measurement precision solves it. You could measure
> perfectly and still get the answer backwards.
>
> Why? Because the standard design compares a tumour to nearby healthy tissue
> from the same patient. Both samples are taken at the same moment, after the
> cancer already exists. So if you find a bacterium is different, there are at
> least four explanations:
>
> One, the bacterium caused the cancer. Two, the cancer changed the tissue and
> the bacterium responded. Three, something else caused both, diet, smoking,
> inflammation. Four, it's an artefact of how the samples were handled.
>
> A single time point cannot distinguish these. Not because the data is bad.
> Because the *design* doesn't contain the information."

## 3. Building data with nothing in it, 2 minutes

> "So my first question wasn't about bacteria at all. It was: if the answer
> were nothing, would our methods say nothing?
>
> This is a strange question to ask, so let me explain why it's the right one.
>
> In a lab, when you run a machine, you run a blank. Pure water through the
> instrument. If the machine reports a result on water, the machine is broken,
> and you'd never trust its readings on real samples. That's a negative
> control.
>
> Statistical methods deserve the same treatment, and they rarely get it.
>
> So I built the blank. Picture a spreadsheet: rows are patients, columns are
> bacteria, each cell is how much of that bacterium was found. Off to the side,
> a few extra columns: what cancer this was, whether the sample is tumour or
> normal, and which lab processed it.
>
> I made two versions.
>
> Version A has real signal. Six cancer types, 360 patients, 150 bacteria. For
> certain bacteria in certain cancers, I nudged the numbers up. I know exactly
> where the signal is and how strong, because I put it there.
>
> Version B has nothing. Same size, same shape, generated by the same code. But
> the bacterial numbers are generated completely independently of the cancer
> labels. There is nothing to find. Not 'probably nothing'. Nothing, by
> construction.
>
> And then the sneaky part. In version B I made which of three labs processed a
> sample correlate with the cancer label, about three quarters of the time.
> That's realistic: hospitals send samples in batches, different centres see
> different patients.
>
> Now there's a shortcut. A pattern-finding algorithm doesn't need biology. It
> can learn 'samples from lab two are usually cancer' and score well. It looks
> like it found something. It found the paperwork."

## 4. What the standard tests do, and why they fail, 2 minutes

> "Two tests dominate this literature. I didn't choose them; I went and looked
> at what gets reported, and reproduced that. If I'd picked my own favourite
> tests and shown they fail, nobody should care.
>
> **Test one, the no-information rate.** If you're predicting six cancer types
> and one is 30 percent of your samples, you can score 30 percent by always
> guessing that one. That's the floor. The test asks: does the model beat the
> floor?
>
> **Test two, label permutation.** Take your spreadsheet and shuffle the
> diagnosis column, like shuffling a deck. Any real link between bacteria and
> diagnosis is destroyed. Run the whole analysis on the shuffled version, many
> times. If the real version scores better than nearly all the shuffled ones,
> you call it significant.
>
> That second one sounds airtight. It's the more sophisticated of the two and
> the one people trust.
>
> Both passed on the dataset containing nothing. Ten seeds out of ten.
>
> And they fail for the *same* reason, which is worth understanding. Neither
> test asks *what* the model learned. Only *whether* it learned.
>
> The model can learn 'lab two means cancer', which beats the no-information
> rate comfortably. And permutation doesn't catch it, because when you shuffle
> the diagnosis column you break the link between bacteria and diagnosis, but
> you *also* break the link between lab and diagnosis. You destroyed the
> shortcut along with the signal. So shuffled performance collapses, the real
> version looks great by comparison, and the test declares victory.
>
> Both tests ask 'is this better than random?' The honest question is 'is this
> better than the dumbest non-biological explanation available?' Those are not
> the same question, and the gap between them is where a lot of this
> literature has been living."

## 5. What does work, 60 seconds

> "Two fixes, neither exotic.
>
> **A confounder baseline.** Instead of comparing your model against random
> guessing, compare it against a model built from the confounder alone. Feed a
> model nothing but the lab identifier, no bacteria at all. Whatever it scores
> is the real floor. Now the question is: do the bacteria add anything on top
> of the paperwork? On the empty data, they don't. The margin goes to zero, and
> the test correctly says no.
>
> **Within-batch validation.** Restructure so the model trains on some labs and
> is tested on a lab it has never seen. Now 'lab two means cancer' is worthless,
> because there is no lab two in the test set. The shortcut is unavailable by
> construction. Real biology should transfer across labs. Paperwork shouldn't.
> On the empty data, performance drops to chance. Correctly."

## 6. The natural control, 2 minutes — the core of the talk

> "Simulation proves the method *can* fail. It doesn't prove it *did* fail on
> real patients. For that I needed an organism where the truth is known
> independently of any microbiome study.
>
> *Helicobacter pylori* is that organism. It's classified as a Group 1
> carcinogen, the same category as tobacco and asbestos. It causes gastric
> adenocarcinoma. Eradicating it reduces cancer incidence. This is about as
> settled as causation gets in cancer biology.
>
> So the test is simple. Take a careful decontaminated dataset, look at stomach
> tumours, and ask: is *H. pylori* more abundant in the tumour than in the
> patient's own adjacent healthy tissue?
>
> If the standard design works, it should be enriched. It's the cause.
>
> It's about eight-fold depleted. Thirty-nine matched pairs. p equals 4.5 times
> ten to the minus four.
>
> Now, this is not a surprising result to a gastroenterologist, and that's the
> point. The mechanism is well described. *H. pylori* colonises the acidic
> gastric mucosa. Chronic infection drives inflammation, then atrophy, then
> intestinal metaplasia, where the stomach lining is replaced with something
> more intestine-like. That new tissue is less acidic and no longer hospitable.
> The organism is displaced by the very process it started.
>
> Cause first, then disappearance. And the cross-sectional snapshot catches only
> the second half.
>
> So for the one organism whose causal role we independently know, the standard
> design recovers the *opposite* of the truth. Not a weaker version. The
> opposite.
>
> That's the finding. Everything else on the poster is setup or consequence."

## 7. Two further ways clean data misleads, 45 seconds

> "Two smaller results, both practical.
>
> **The detection floor.** I measured how big an effect had to be before my
> method could reliably find it. About 0.8 log units at this sample size. Below
> that, real effects are missed most of the time. Which means when I report a
> null, the honest statement is 'no signal above 0.8 log units', not 'no
> signal'. Most nulls in this literature don't come with that number, which
> makes them uninterpretable.
>
> **Taxonomic redundancy.** Bacteria are classified in a nested hierarchy:
> phylum contains class contains order contains family contains genus. If a
> signal exists at one level, it often appears at all of them. Whether you count
> those as one finding or five is a choice, and it's rarely stated. On identical
> data at an identical significance threshold, my discovery count varied 5.7-fold
> depending on that choice alone."

## 8. Germline anchoring, 2 minutes

> "So cross-sectional abundance can't establish direction. What can?
>
> You need a variable that the disease cannot possibly have affected. Genes
> qualify. You inherited them at conception, decades before any tumour. A tumour
> cannot reach back and change your genome.
>
> And it turns out some genetic variants influence which bacteria you host. Not
> strongly, but measurably. So nature has effectively run a randomised trial:
> at conception, people were randomly assigned slightly more or slightly less
> of a given bacterium.
>
> Compare cancer rates between those groups, and the direction is guaranteed by
> chronology. Genes came first. This is called Mendelian randomisation.
>
> I ran it for 211 gut bacteria against colorectal cancer, using a Finnish
> cohort of about 11,800 cases. Then I repeated it in a completely independent
> cohort of 100,204 cases.
>
> The result is a null, and I want to be precise about what kind.
>
> Nothing survived correction for multiple testing in the large cohort. Every
> lead from the smaller cohort failed to replicate. The single hit that did
> survive correction in the first cohort was *Cyanobacteria*, and it doesn't
> hold up: its pleiotropy test sits right at p 0.05, one of its genetic
> instruments sits at NOS2, a gene directly implicated in colorectal
> carcinogenesis, which violates the core assumption that the gene affects
> cancer *only* through the bacterium, and its direction flipped in the
> replication cohort.
>
> There's also a nice validation. A *Bifidobacterium* signal looked promising
> until I removed a single genetic variant, rs182549. That variant is at the
> lactase gene, the one that determines whether you can digest dairy as an
> adult. It's a notorious trap in this field, because dairy intake affects both
> gut bacteria and colorectal cancer risk independently. Removing it collapsed
> the signal, odds ratio from 1.21 to 1.08, p from 0.02 to 0.27. The pipeline
> caught it without being told to look."

## 9. The residual signal, 60 seconds

> "But the null isn't the whole story.
>
> I compared the two independent cohorts, not on which bacteria were
> significant, but on which *direction* each pointed. Protective or harmful.
>
> Across 210 bacteria tested in both, they agreed on direction 61.4 percent of
> the time. Under pure noise you'd expect 50. That's a binomial p of about
> 0.001.
>
> So there is real causal signal from the gut microbiome on colorectal cancer.
> It's just diffuse. Spread thinly across many organisms, none individually
> strong enough to survive multiple testing at current sample sizes.
>
> That matters for how you'd act on it. If the causal architecture is one
> organism, you develop a targeted therapy. If it's a hundred organisms each
> contributing a little, single-strain products are the wrong bet and
> community-level interventions are the right one. That's a very different
> investment thesis, and the evidence currently points at the second."

## 10. Close, 45 seconds

> "Let me land three things.
>
> First, decontamination is necessary and not sufficient. The field solved
> 'is it there' and the harder question is still open.
>
> Second, for the one organism where we know the truth independently, the
> standard design gives the wrong direction. That's not a hypothetical failure
> mode. It's a demonstrated one, in a careful dataset.
>
> Third, and this generalises well past microbiome work: a test that cannot
> fail is not a test. If you never run your method on data you know contains
> nothing, you have no idea what its output means when it says something.
>
> That last one isn't a microbiology point. It's a machine learning point, and
> it applies to imaging models, risk scores, and most omics.
>
> Happy to take anything."

---

# VERSION 4 — Q&A bank

## The ones you will actually get

**"So which bacterium causes colorectal cancer?"**
> "None that I can name, and I think that's the finding. No single taxon
> survives correction in a hundred thousand cases. What survives is a diffuse
> directional signal, 61 percent agreement across independent cohorts against
> 50 under noise. That's consistent with many small contributions rather than
> one driver."

**"Isn't the H. pylori result just well-known gastric pathology?"**
> "Yes, entirely. That's why it works as a control. I'm not claiming a
> discovery about stomach cancer. I'm using an organism whose causal role is
> settled to test whether the *method* recovers it. It doesn't. If the result
> were novel it would be a worse control, because then we'd be arguing about
> the biology instead of the method."

**"Aren't you just saying everyone's findings are wrong?"**
> "No, and I'd resist that framing. I'm saying the tests as reported cannot
> distinguish real signal from artefact. Some of that work may be entirely
> correct. We can't tell from what's on the page. That's a weaker and more
> defensible claim, and it's the one the data supports."

**"Your simulation is artificial. Real data isn't like that."**
> "Agreed, and that's why the simulation isn't the argument. It establishes
> that the failure mode is possible and that the standard tests can't see it.
> The real-data half is the H. pylori inversion, in a decontaminated published
> cohort that passed every control I ran. The simulation tells you what to look
> for; the real data tells you it's there."

**"Why random forests? Wouldn't a better model fix this?"**
> "It would make it worse. A more flexible model finds the shortcut more
> efficiently. The problem isn't model capacity, it's that the validation
> doesn't ask what was learned. I used a deliberately unremarkable model
> precisely so the finding can't be about the model."

**"Only 39 pairs. Isn't that small?"**
> "It is, and it's the main limitation of that panel. It's also every matched
> gastric pair in the dataset. The effect is large, about eight-fold, and
> consistent in direction across pairs. I'd very much like an independent
> matched-pair gastric cohort, and I haven't found one. If you know of one, I'd
> like to hear about it after this."

**"Isn't 61 percent quite weak?"**
> "It is weak, and I present it as weak. It is not, however, ambiguous:
> binomial p of 0.001 over 210 taxa, and it holds under every rule I tried for
> handling ties. The claim is only that the signal is non-zero and diffuse, not
> that it's actionable."

**"Mendelian randomisation has strong assumptions."**
> "It does, and I'd rather state them than defend them. The instruments are
> weak, from a p<1e-5 threshold rather than genome-wide significance.
> Instrument strength is the binding constraint, not outcome power. I've done
> no colocalisation and no reverse-direction MR. And genus-level instruments
> may not correspond to biologically meaningful units. All of that is in the
> limitations. My view is that species-level instrumentation is the only thing
> that would meaningfully change this null."

## The hostile ones

**"This is a negative result. Why is it a poster?"**
> "Because the negative result isn't the contribution. The contribution is that
> the standard tests pass on data containing nothing, ten out of ten, and that
> a known carcinogen appears on the wrong side of the comparison. Those are
> positive findings about the methods. The null is what falls out afterwards."

**"Who peer-reviewed this?"**
> "Nobody yet, it's under preparation. The code, the committed results, the
> figure generation and an append-only corrections log are all public, and
> continuous integration re-derives every number in the manuscript on each
> commit. If a value ever stops matching, the badge goes red. That's not peer
> review, but it's more than most posters here can offer, and I'd rather you
> check than take my word."

**"Have you found errors in your own work?"**
> "Yes, twelve recorded corrections, all public. One mattered: an effect size
> in my notes, odds ratio 0.83 at p 0.026, did not exist in my saved data. The
> real value was 0.95 at p 0.58, a null. Root cause was an analysis that ran
> interactively and was never saved, plus a results table that stored a
> formatted display string instead of the underlying number. Both are fixed and
> both are written up. Given what this poster argues, hiding that would be
> indefensible."

---

# THINGS TO NEVER SAY

| Do not say | Say instead |
|---|---|
| "Alistipes is protective" | "Alistipes points protective in three of three tests and is nominally significant in two. The colorectal FinnGen test is null at p 0.578. It's a lead, not a finding" |
| "The field is wrong" | "The tests as reported cannot distinguish" |
| "We found no effect" | "No signal above a detection floor of 0.8 log units" |
| "H. pylori doesn't cause gastric cancer" | "H. pylori causes gastric cancer and is nonetheless depleted at the tumour site" |
| "61 percent, so the microbiome matters" | "61.4 percent against 50 under noise, robust to the tie rule, indicating a real but diffuse signal" |
| "2.5x chance from confounding" (unqualified) | "2.5 times chance as confounding rose to 0.95" |
| Any p-value on the three-test triangulation | "82 of 210 taxa agree in direction across all three tests, against about 52 expected. Descriptive only, the tests aren't independent" |

---

# NUMBERS CARD — keep in your pocket

| Quantity | Value |
|---|---|
| Standard tests passing on zero-signal data | 10 of 10 seeds |
| Confounding-only accuracy | 2.46x chance at 0.95 confounding; 1.94x at 0.8 |
| H. pylori depletion | ~8-fold; 39 matched pairs; difference −0.994 log units; p = 4.5e−4 |
| Within-tissue AUC | COAD 0.50, ESCA 0.58, STAD 0.66, HNSC 0.69 |
| Detection floor | ~0.8 log units at n = 360 |
| Nesting rule effect | up to 5.7-fold on discovery count |
| MR screen | 211 taxa; FinnGen 11,790 cases; replication 100,204 cases |
| FDR survivors in the large cohort | zero |
| Cross-cohort directional agreement | 61.4%, 129/210, binomial p = 0.0011 |
| Alistipes | protective direction 3/3; nominal 2/3; FinnGen colorectal p = 0.578 |
| Cyanobacteria | Egger intercept p = 0.05; NOS2 instrument rs2314810; direction flipped on replication |
| Bifidobacterium / LCT | OR 1.21 → 1.08, p 0.02 → 0.27 on removing rs182549 |
| Cohort | TCMA, 611 samples, 14,492 taxa, 5 TCGA projects, 3 centres |
