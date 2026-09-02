---
title: "Instrumental Variables Before and After an Event"
subtitle: "Identification with Endogenous Exposure"
author: "Sriteja Burla"
date: "August 2026"
lang: en-US
---


## Summary

Suppose a hospital is acquired. Did the acquisition change the effect of receiving care there? A direct comparison with patients treated elsewhere is difficult because hospital destination is not random. Ambulance dispatch offers a possible instrument: otherwise similar patients may be assigned to companies whose hospital preferences shift which hospital initially receives them (Doyle et al. 2015).

The tempting calculation is to estimate the ambulance IV before and after acquisition and subtract the Wald ratios. It is not yet an event effect. Each ratio can identify a local effect for patients whose destination responds to the instrument in that period, but their difference requires the same target population, comparable complier populations, and a counterfactual path for the LATE absent acquisition. Exclusion must hold in the observed periods and on that unobserved path too.

Under the benchmark used here, the same population complies before the event and under both post-event states, and the counterfactual post-event LATE equals the comparable pre-event LATE. The event effect then equals the standardized post-event Wald ratio minus the standardized pre-event Wald ratio. Sensitivity analysis can examine less restrictive assumptions about the complier populations or the missing counterfactual; throughout, I separate what the period-specific IVs identify from the assumptions that turn their difference into an event effect.

## 1. The empirical setting

### 1.1. Roles of the event and instrument

I use a hospital acquired at a known date as the running example: the acquisition is the event, and emergency patients are the units whose hospital exposure may change. Let $D=1$ denote initial care at the focal hospital, while $D=0$ denotes a prespecified alternative exposure—for example, initial care at one comparison hospital or at hospitals drawn from a fixed eligible-hospital distribution. Let $Z\in\{0,1\}$ denote a prespecified encouragement derived from ambulance-company assignment or company preferences.

Ambulance assignment and acquisition do different work: the former shifts hospital destination within a period, whereas the latter may change care at both the focal hospital and the alternatives against which that care is compared. Thus, $Z$ identifies exposure effects within each period; reading their change as an event effect calls for further assumptions across periods.

Fuzzy or instrumented DiD gets its instrument from a policy or group change (de Chaisemartin and D'Haultfœuille 2018; Miyaji 2024). In this design, by contrast, a within-period IV identifies exposure effects, while the event changes the environment across which those effects are compared.

\Needspace{0.62\textheight}

**Figure 1. The event and the instrument play different causal roles**

![](figures/figure1_nested_event_roles.pdf){width=88%}

*Notes:* The IV shifts exposure within a period. The event may change both the exposure process and the effect of exposure. Exclusion rules out an effect of $Z$ on $Y$ through channels outside the prespecified exposure contrast.

Hospitals are only one example. The same problem arises when a workplace policy prompts workers to sort across establishments or when families choose among schools that adopt programs. An institutional reform that changes case assignment across courts creates it too. In each setting, an instrument addresses exposure choice within a period, and an event changes the environment in which that exposure operates.

### 1.2. Target population and intervention definitions

Define the **target population** by a rule that remains fixed across the event and throughout the analyzed periods. In the hospital application, it might include all eligible emergency patients originating in prespecified pickup areas around the hospital, whether they receive care at the focal hospital or are ultimately taken elsewhere.

Membership must be determined before ambulance assignment; otherwise, realized destination—the exposure moved by the instrument—defines the study population. Predetermined market, cohort, or source-population indicators can refine the target, but realized-provider fixed effects condition on an endogenous response and may absorb the first stage or select patients on unobserved determinants of outcomes. Konetzka, Yang, and Werner (2019) discuss the related problem of using a patient-level instrument to study a provider-level attribute.

The interventions must retain the same meaning across periods: the two values of $Z$ should represent the same encouragement, while $D=1$ and $D=0$ should continue to describe the same exposure contrast. Pooling several hospitals in the comparison group makes this especially difficult because a change in their mixture may alter the meaning of $D=0$. Fixing the alternative-provider distribution addresses that problem. If this is infeasible, the analysis must defend version irrelevance or define a richer exposure.

The data need not form a panel. If eligibility or outcome observation changes over time, the records may instead be repeated cross sections. In either case, state the observation rule alongside the target population. Any exposure-induced survival, attrition, or missingness changes which units contribute an observed outcome and requires a separate selection analysis.

## 2. Identification of the event effect

### 2.1. Period-specific local average treatment effects

Let $t=0$ denote a pre-event period and $t=1$ a post-event period. For predetermined covariates $X=x$, the conditional Wald ratio is

\[
\beta_t(x)=
\frac{E[Y\mid Z=1,X=x,t]-E[Y\mid Z=0,X=x,t]}
{E[D\mid Z=1,X=x,t]-E[D\mid Z=0,X=x,t]}.
\tag{1}
\]

Write the numerator as $\rho_t(x)$ and the denominator as $\pi_t(x)$. Under the usual binary-IV conditions—independence, exclusion, monotonicity, positive probability for both instrument values, and relevance—$\beta_t(x)$ is the exposure LATE for period-$t$ compliers whose outcomes are observed (Imbens and Angrist 1994). Define $Z=1$ as the value intended to increase exposure and keep that coding in every period.

These IV conditions must hold in the analyzed sample for every period used in the event design; when outcome observation is selective, or the estimand is meant to cover all compliers in the target population, the additional conditions in Section 2.2 are also needed. Report the instrument propensity, reduced form, first stage, and strength evidence separately by period because a pooled first stage can conceal a weak component.

### 2.2. Standardization and outcome observation

Raw Wald ratios can change even when the conditional LATEs do not, simply because the observed covariate distribution of compliers changes. Choose a reference distribution $F_X^\star$ before examining outcomes, then average the conditional ratios over the part of its domain with common support:

\[
\beta_t^\star=E_{F_X^\star}[\beta_t(X)],
\qquad
\Delta^{IV,\star}=\beta_1^\star-\beta_0^\star.
\tag{2}
\]

The covariate distribution of pre-event compliers is a useful default. Given the observed pre-event distribution $F_{X0}$, the corresponding reference measure is proportional to $\pi_0(x)dF_{X0}(x)$ and should be restricted to strata that support both period-specific ratios. Another policy-relevant distribution is equally coherent, provided researchers report that choice and any trimming.

Standardization puts the pre- and post-event Wald ratios on the same observed covariate distribution, but does not, by itself, extend the estimand beyond compliers whose outcomes are observed. That extension requires the conditional mean exposure effects of target-population compliers to equal those of observed compliers, as well as positive probability of outcome observation in every target stratum; exclusion must also hold for the target population, not only for the analyzed records. Selective survival or attrition will often require an explicit model or sensitivity analysis.

### 2.3. Complier populations and the post-event counterfactual

Standardization fixes the observed covariate distribution, but it doesn't guarantee that $\Delta^{IV,\star}$ compares the same compliers. The benchmark assumes that the same units would comply before the event and after it under both the observed and counterfactual states; a weaker approach allows the complier population to change while requiring equality of the relevant conditional mean exposure effects across those populations. Sections B.2 and B.3 state the two sets of conditions.

Even with comparable complier populations, a second question remains: what would the post-event LATE have been at the same date without the event? Let $L_1^{E,\cap,\star}$ denote the standardized post-event LATE under the event and $L_1^{0,\cap,\star}$ that missing counterfactual; both refer to units who comply in all three relevant states and use the reference covariate distribution. The event estimand is

\[
\theta^{\cap,\star}=L_1^{E,\cap,\star}-L_1^{0,\cap,\star}.
\tag{3}
\]

The data reveal only the first term after the event. The simplest assumption equates the missing term with the comparable pre-event LATE; with several pre-event estimates, researchers can instead extrapolate the pre-event LATE path linearly or bound plausible departures from it. This assumption concerns the path the LATE would have followed, not a trend in average untreated outcomes.

Exclusion is needed in the counterfactual post-event state too: holding exposure fixed, $Z$ has no effect on the outcome that would be observed without the event. Since that state is unobserved, the case for this assumption must come from the institutional design. No anticipation, meanwhile, gives the pre-event estimates their counterfactual meaning; announcements, preparation, or early changes in ambulance behavior may contaminate periods close to the acquisition date, in which case use an earlier reference period or redefine the transition window.

\Needspace{0.58\textheight}

**Figure 2. The observed IV path leaves one counterfactual LATE missing**

![](figures/figure2_missing_late_path.pdf){width=88%}

*Notes:* The observed pre-event and post-event Wald ratios identify period-specific local effects under their respective IV conditions. An event interpretation also requires a comparable complier population and a specification for how the post-event LATE would have evolved without the event. The dashed line is schematic: applications may hold the pre-event LATE fixed, extrapolate its earlier path, or bound plausible counterfactual values.

Given the period-specific IV conditions and the event assumptions above, the guide's core result is

\[
\boxed{\theta^{\cap,\star}=\beta_1^\star-\beta_0^\star.}
\tag{4}
\]

The observed data provide the two Wald ratios in Equation (4). Everything that turns their difference into the event effect comes from the restrictions on outcome observation and complier membership and from the assumed counterfactual path of the LATE.

### 2.4. Identified causal objects

Table 1 separates the claims supported at each stage of the design.

\Needspace{0.38\textheight}

**Table 1. What the design supports at each stage**

\begingroup
\renewcommand{\arraystretch}{1.30}

| Maintained conditions | Supported interpretation |
|---|---|
| Outcome and exposure means by instrument value are well defined | Reduced forms and first stages |
| Period-specific IV conditions hold | A LATE of exposure for each period's observed compliers |
| Common support and standardization hold, and conditional mean effects are equal for observed and target-population compliers | Average conditional exposure effects for target-population compliers, standardized to the chosen reference distribution |
| The same complier population is used across periods and event states; no anticipation and exclusion hold; and the counterfactual post-event LATE is correctly specified | Change in the standardized LATE caused by the event for that complier population |

\endgroup

The final row is demanding. If it is hard to defend, the standardized change in period-specific LATEs still has a narrower but useful interpretation: exposure effects changed for the complier populations reached by the instrument. Researchers should settle on that interpretation before estimation and state what they will report if the stronger assumptions prove unconvincing.

## 3. Decomposition and sensitivity analysis

### 3.1. Diagnostic decomposition

A standardized IV change need not be the event effect: the LATE could have moved even without the event, the instrument may reach a different complier population, or the observed records may cover a different part of the target population. The decomposition separates these possibilities:

\[
\Delta^{IV,\star}
=\theta^{\cap,\star}
+\Delta^{0,\star}
+B^{C,\star}
+B^{O,\star}.
\tag{5}
\]

Here $\Delta^{0,\star}$ is the change in the standardized LATE that would have occurred without the event for the common complier population. $B^{C,\star}$ captures changes in that population, while $B^{O,\star}$ captures changes in the difference between observed compliers and compliers in the target population; together, the terms lay bare what must be assumed before a change in period-specific IV estimates can be called an event effect. Use the decomposition only after defending the period-specific IV conditions and intervention definitions.

Each term calls for different evidence. Historical IV estimates can inform $\Delta^{0,\star}$ if they use the same instrument, exposure, target population, and assumptions about complier comparability; covariate balance and standardization to the same distribution help diagnose observable composition, while first-stage changes reveal changes in complier-group size, although similar first stages can hide substantial entry and exit. Sampling rates, attrition patterns, and audits of outcome observation speak to $B^{O,\star}$. No one diagnostic settles the decomposition.

The three gaps can instead be bounded jointly. If their absolute values are bounded by $h_0$, $h_C$, and $h_O$, then

\[
\theta^{\cap,\star}
\in
\left[\Delta^{IV,\star}-H,\ \Delta^{IV,\star}+H\right],
\qquad H=h_0+h_C+h_O.
\tag{6}
\]

Because the discrepancies may occur together, the sensitivity analysis should vary them jointly; Rambachan and Roth (2023) provide a template for calibrating restrictions on the counterfactual LATE path. Inference must be developed for the resulting IV estimand.

### 3.2. Simulation design and results

I use Figure 3 to make the three problems concrete. The simulation follows four repeated cross sections around one event in which the IV is valid in every period and the event causes a true change in LATE of one in each scenario. Across scenarios, it separately varies observed covariate composition, unobserved complier types, and the counterfactual LATE path without the event.

\Needspace{0.56\textheight}

**Figure 3. Each assumption addresses a different source of bias**

![](figures/figure3_standardization_and_continuation.png)

*Notes:* Markers are means across 2,000 replications; vertical lines mark the 2.5th and 97.5th percentiles. Each period has 2,500 observations. The dashed line gives the true change in LATE caused by the event. “Standardized” averages over the pre-event complier covariate distribution. “Pre-event LATE held fixed” uses the comparable pre-event LATE as the counterfactual post-event value; “Linear extrapolation” uses three pre-event LATEs. The replication materials document the data-generating process, fixed seed, and outputs.\footnote{Replication materials: \url{https://github.com/TejaBurla89/iv_nested_event}.}

In the benchmark case, the same units comply in all three states, and all three estimators average about one. When only the observed distribution of $X$ changes, the unstandardized estimate rises to 1.75; standardizing to the pre-event complier distribution brings it back to 1.00.

That correction goes no further. Turnover in unobserved complier types leaves the standardized estimate near 1.80 even though the observed covariate distribution does not change. In the last scenario, the LATE would have risen without the event: holding the pre-event value fixed gives 1.40, whereas linear extrapolation from three pre-event LATEs gives 0.99. An application therefore needs evidence for each claim, and it will usually come from different sources.

## 4. Multiple periods and staggered events

### 4.1. Event-time local average treatment effects

With more periods, the two Wald ratios become an event-time path of standardized LATEs whose pre-event portion shows how the local effect was evolving and whether the first stage or common support changed.

These patterns can guide the choice of a post-event counterfactual, but they do not identify what would have happened beyond the event. Low-powered pretrend tests provide limited reassurance, and choosing an extrapolation because a particular specification passes such a test changes the subsequent inference (Roth 2022); the design should therefore specify the reference period, transition window, counterfactual assumption, and sensitivity range in advance.

### 4.2. Cohort-specific estimation and aggregation

Staggered timing creates another source of changing composition because different cohorts may contribute at different event times. I would start with cohort-by-event-time estimates so those changes remain visible, following the group-time discipline used in modern staggered designs (Callaway and Sant'Anna 2021). Within a cohort, calculate conditional reduced forms and first stages for every institution or other event unit; form Wald ratios in covariate strata with adequate support, standardize them, and only then average across event units.

The order matters because it preserves the prespecified aggregation weights; pooling reduced forms and first stages before taking the ratio creates first-stage-dependent weights and generally defines a different question. Aggregate the cohort-specific effects using prespecified weights. Comparisons across event time are cleaner when the same cohorts and weights enter every horizon, but each cohort-period estimate still requires adequate first-stage and covariate support, adequate outcome observation, and credible assumptions relating its complier population to the others. If support disappears, report the share of the reference distribution no longer represented, then either defend an extrapolation or relabel the estimand for the supported set.

Stacking duplicates some original observations across cohort comparisons, but those copies carry the same sampling information; retain the original cluster identifiers and rebuild the stacks within a bootstrap. Calendar-time shocks also deserve attention when they change the exposure contrast, compliance behavior, observation, or instrument validity.

## 5. Estimation, inference, and reporting

### 5.1. Estimation

Before looking at outcomes, I would fix the event date and transition window, then define the target population, instrument construction, and exposure contrast. The same plan records the covariates and support rule, the reference distribution, the outcome-observation requirements, and the assumptions used to compare complier populations and construct the counterfactual post-event LATE. A staggered design also fixes the weights within and across cohorts.

Equation (2) then determines the calculation. Estimate the conditional reduced forms and first stages, form Wald ratios wherever both periods have adequate support, and average them using the reference weights. Regression software handles the regressions; it cannot choose the target population or weighting scheme. I would reproduce the result directly, on the same sample and with the same weights, as a check.

When the instrument is generated, construct it from predetermined information that the event cannot affect, fixing the training population and dates, construction rule, sign, cutoff, omitted observations or clusters, and treatment of unsupported inputs. If data from the main analysis help estimate the score, rebuild it with the rest of the estimator during resampling; leaving out the focal observation or cluster prevents mechanical contribution to its own instrument. None of these construction choices substitutes for an institutional argument for independence and exclusion.

### 5.2. A numerical example

Table 2 works through the calculation for two covariate groups whose reference weights are 0.60 and 0.40. The standardized LATE is 1.68 before the event and 2.78 afterward, a change of 1.10; if all the conditions in Table 1 hold—including the assumption that the counterfactual post-event LATE equals the pre-event LATE—then 1.10 is the event estimate.

\Needspace{0.43\textheight}

**Table 2. From means by instrument value to the standardized contrast**

| Period | Group | Reference weight | Reduced form | First stage | Wald ratio |
|---|---:|---:|---:|---:|---:|
| Pre | A | 0.60 | 0.36 | 0.20 | 1.80 |
| Pre | B | 0.40 | 0.24 | 0.16 | 1.50 |
| Post | A | 0.60 | 0.58 | 0.20 | 2.90 |
| Post | B | 0.40 | 0.40 | 0.16 | 2.50 |

The pre-event standardized LATE is $0.60(1.80)+0.40(1.50)=1.68$; after the event, it is $0.60(2.90)+0.40(2.50)=2.78$. Table 1 governs what the 1.10 difference means. Under its first three rows, the difference is the standardized change in period-specific LATEs.

### 5.3. Inference

The components of an event estimate are generally dependent because they may share observations, reference weights, counterfactual estimates, cohort stacks, or a generated instrument; a ratio also becomes unstable when its first stage is weak. Under strong identification, I would therefore conduct inference on the final event contrast, using, for example, a cluster bootstrap that resamples the original independent units and re-estimates every component, including the reference distribution, cohort stacks, and generated scores.

Weak first stages call for identification-robust inference on the final event estimand. Section E.3 gives an executable conservative projection in one deliberately narrow case: one pre-event period, one post-event period, fixed finite strata, fixed external weights, and a counterfactual post-event LATE equal to the pre-event LATE. Simultaneous intervals for the reduced forms and first stages allow the final set to become unbounded when the data cannot exclude a zero denominator. Estimated weights, selected support, and staggered aggregation fall outside that construction and need inference tailored to the larger estimator; if none is available, the honest fallback is to report the reduced-form and first-stage moments and label conventional intervals as strong-identification approximations.

Choose clusters to match the source of assignment and dependence—often the ambulance company, market, hospital, or some combination—rather than merely the patient row. Few-cluster and weak-IV problems are distinct and require separate remedies. When its conditions hold, a wild cluster bootstrap may improve reduced-form and first-stage inference with few clusters (Cameron, Gelbach, and Miller 2008), but resampling a conventional ratio statistic does not make that statistic identification robust.

### 5.4. Reporting

A reader should be able to reconstruct each event estimate. Start with the outcome and exposure means by instrument value, then give the reduced forms, first stages, reference weights, and standardized period-specific ratios. Report support losses and observation rates separately, including the share of the reference distribution retained. The account of the counterfactual LATE should identify its inputs; for a staggered design, give the cohort weights and clustering structure. If the instrument is generated, describe the training sample and say whether inference rebuilds its construction.

Table 1 should also govern the label attached to the final estimate. I would reserve “change in the standardized LATE caused by the event for units who would comply in all three states” for cases supporting the full identification argument; with weaker evidence on complier comparability or the counterfactual LATE, use “standardized change in period-specific LATEs.” Statistical precision indicates how much the sample reveals about the estimand actually chosen.

## 6. Conclusion

Two separately valid IV estimates do not, by themselves, identify an event effect.

The within-period instrument deals with selection into exposure. It does not choose a common target population. Nor can it establish that the complier populations are comparable or reveal the LATE that would have prevailed after the event in its absence. Those parts of the argument come from assumptions. If they are credible, the standardized difference in Wald ratios identifies the change in the standardized LATE caused by the event. If they are not, the period-specific LATEs and sensitivity ranges are the more defensible result—they preserve what the IV estimates reveal without hiding what an event interpretation would require.

\vspace{2.2em}

\noindent\rule{0.34\linewidth}{0.4pt}

\vspace{1.0em}

**Suggested citation**

Burla, Sriteja. 2026. *Instrumental Variables Before and After an Event: Identification with Endogenous Exposure*. Version 1.0, August 2026.

\newpage

# Technical Supplement

## A. Setup and within-period IV

Fix an institution, or some other unit at which the event occurs, and let $\mathcal U$ be its target population, with membership in $\mathcal U$ determined before the instrument and event. The index $a\in\{0,E\}$ records whether the event has occurred by the post-event date: at $t=1$, $a=E$ is observed and $a=0$ is counterfactual. Under no anticipation (Section B.2), the two paths coincide at $t=0$, so the observed states are $a_0=0$ and $a_1=E$.

For individual $i$ in period $t$, let $D_{it}(a,z)\in\{0,1\}$ denote potential exposure under binary instrument value $z$ and $Y_{it}(a,d)$ the potential outcome under exposure $d$. The notation builds in exclusion: after fixing exposure, the instrument has no further effect on the outcome. In the observed periods, exclusion belongs to the period-specific IV argument, while extending it to the counterfactual post-event state without the event takes a separate assumption.

Keep the instrument coding fixed across periods, with $z=1$ denoting the value intended to increase exposure; monotonicity then requires $D_{it}(a,1)\geq D_{it}(a,0)$. Define the complier set and exposure effect by

\[
\mathcal C_t(a)=\{i\in\mathcal U:D_{it}(a,1)>D_{it}(a,0)\},
\qquad
\tau_{it}(a)=Y_{it}(a,1)-Y_{it}(a,0).
\tag{S.1}
\]

Under the usual binary-IV conditions, Equation (1) identifies the conditional mean exposure effect among that period's observed compliers—the standard LATE result (Imbens and Angrist 1994); Section 2.1 discusses validity and reporting period by period.

The indicator $O_{it}$ records whether the outcome of an eligible unit is observed. Conditional on predetermined covariates, the benchmark assumes that neither the instrument nor exposure changes whether that outcome is observed; if survival, attrition, or missingness responds to either one, a separate selection argument is needed.

## B. Standardization and event identification

### B.1. Standardization and outcome observation

Let $\mathcal X^\star$ contain the covariate values supported in both observed periods, with reference distribution $F_X^\star$ as defined in Section 2.2.

Moving from observed compliers to the target population requires their conditional mean effect to equal the corresponding mean among all compliers in that population:

\[
E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i=x,O_{it}=1]
=E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i=x].
\tag{S.2}
\]

When Equation (S.2) holds, the standardized Wald ratio identifies the average conditional exposure effect for target-population compliers, standardized to $F_X^\star$:

\[
\beta_t^\star
=E_{F_X^\star}\!\left[
E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i]
\right].
\tag{S.3}
\]

Equation (S.2) restricts the mean alone; transporting the joint distribution of potential exposures and outcomes would take a stronger assumption. Either version also requires positive observation probability wherever $F_X^\star$ assigns weight.

### B.2. A common complier population and the post-event counterfactual

Let $\mathcal C^\cap$ be the intersection of the complier sets before the event and under the two possible post-event states:

\[
\mathcal C^\cap
=\mathcal C_0(0)\cap\mathcal C_1(0)\cap\mathcal C_1(E).
\tag{S.4}
\]

The main result assumes that all three sets equal $\mathcal C^\cap$ within every value of $X$ receiving positive reference weight and defines, for this population,

\[
L_t^{a,\cap,\star}
=E_{F_X^\star}\!\left[
E[\tau_{it}(a)\mid i\in\mathcal C^\cap,X_i]
\right].
\tag{S.5}
\]

No anticipation means that potential exposure and outcomes before the event cannot depend on the future event path:

\[
D_{i0}(E,z)=D_{i0}(0,z),
\qquad
Y_{i0}(E,d)=Y_{i0}(0,d)
\]

for every relevant $z$ and $d$. These equalities justify no-event notation for the observed pre-event quantities in Equations (S.4)--(S.6), but they do not supply the missing $L_1^{0,\cap,\star}$; a prespecified rule must map the comparable pre-event LATE history into it. The simplest rule assumes $L_1^{0,\cap,\star}=L_0^{0,\cap,\star}$.

Exclusion must also hold in the counterfactual post-event state without the event, an unobserved state for which the restriction has to be assessed from institutional evidence.

**Identification result.** Suppose period-specific IV validity and Equation (S.2) hold. If the three complier sets are equal, there is no anticipation, exclusion extends to the counterfactual post-event state, and the counterfactual post-event LATE equals the pre-event LATE, then

\[
\beta_1^\star-\beta_0^\star
=L_1^{E,\cap,\star}-L_0^{0,\cap,\star}
=L_1^{E,\cap,\star}-L_1^{0,\cap,\star}
=\theta^{\cap,\star}.
\tag{S.6}
\]

This equality completes the identification argument; under a different assumption about the counterfactual LATE, $\beta_0^\star$ is replaced by a prespecified value constructed from the pre-event LATE history.

### B.3. Allowing the complier population to change

The three complier sets need not contain the same units because the algebra requires only equality of the relevant mean effects. To state that weaker condition precisely, write, for any group $G$,

\[
M_{t,a}^{G}(x)
=E[\tau_{it}(a)\mid i\in G,X_i=x].
\]

The two-period design requires the following conditional mean equalities:

\[
\begin{aligned}
M_{0,0}^{\mathcal C_0(0)}(x)
&=M_{0,0}^{\mathcal C_1(0)}(x),\\
M_{1,0}^{\mathcal C_1(0)}(x)
&=M_{1,0}^{\mathcal C^\cap}(x),\\
M_{1,E}^{\mathcal C_1(E)}(x)
&=M_{1,E}^{\mathcal C^\cap}(x).
\end{aligned}
\tag{S.7}
\]

The first equality assigns the pre-event mean effect to the population that would comply after the event if the event did not occur, and the assumption about the counterfactual LATE then carries that population's effect forward. If its counterfactual post-event mean equals its pre-event value, $M_{1,0}^{\mathcal C_1(0)}(x)=M_{0,0}^{\mathcal C_1(0)}(x)$.

The second equality assigns the counterfactual post-event mean effect to $\mathcal C^\cap$, while the third does the same for the observed post-event mean; neither requires the complier populations to contain the same units.

All three restrictions apply wherever $F_X^\star$ assigns weight, and $\mathcal C^\cap$ must have positive probability there. When group membership differs, equal conditional mean effects need a substantive justification.

## C. Decomposition and sensitivity

### C.1. Decomposition

Let $m_t^O$ denote the standardized average conditional effect among observed compliers in period $t$; under period-specific IV validity, $m_t^O=\beta_t^\star$. Let $m_t^U$ be the corresponding standardized average among all period-$t$ compliers in the target population. The three discrepancies in Equation (5) are

\[
\begin{aligned}
\Delta^{0,\star}
&=L_1^{0,\cap,\star}-L_0^{0,\cap,\star},\\
B^{C,\star}
&=(m_1^U-L_1^{E,\cap,\star})-(m_0^U-L_0^{0,\cap,\star}),\\
B^{O,\star}
&=(m_1^O-m_1^U)-(m_0^O-m_0^U).
\end{aligned}
\tag{S.8}
\]

Adding and subtracting these intermediate effects gives Equation (5): the assumption that the LATE would have remained unchanged without the event sets $\Delta^{0,\star}=0$; equality of the complier sets—or the conditional mean equalities in Section B.3—sets $B^{C,\star}=0$; and Equation (S.2) sets $B^{O,\star}=0$.

### C.2. Sensitivity region

Let $\mathcal H$ collect the prespecified plausible triples $(\delta,b_C,b_O)$ for these discrepancies; the corresponding set for the event effect is

\[
\Theta(\mathcal H)
=\left\{
\Delta^{IV,\star}-\delta-b_C-b_O:
(\delta,b_C,b_O)\in\mathcal H
\right\}.
\tag{S.9}
\]

Equation (6) is the symmetric rectangular case, while $\mathcal H$ may more generally impose asymmetric bounds, signs, or dependence across the discrepancies. Historical IV estimates can inform the counterfactual LATE component when they use the same instrument, exposure, target population, reference distribution, and assumptions about complier comparability.

## D. Multiple periods and staggered events

With multiple periods, the two-period argument applies separately to every cohort and event-time component; what changes is the aggregation.

Let $c$ index event cohorts, $h$ institutions or other event units, and $k$ event time. Write a component's conditional reduced form and first stage as $\rho_{hck}(x)$ and $\pi_{hck}(x)$; given the reference distribution $F_{X,hc}^\star$, its observed period-specific standardized LATE is

\[
\beta_{hck}^\star
=E_{F_{X,hc}^\star}\!\left[
\frac{\rho_{hck}(X)}{\pi_{hck}(X)}
\right].
\tag{S.10}
\]

For a reported post-event horizon, define $\mathcal C_{hck}^\cap=\mathcal C_{hck}(0)\cap\mathcal C_{hck}(E)$. The prespecified counterfactual LATE without the event, $g_{hck}$, is constructed from comparable pre-event LATEs for the same complier population and reference distribution. Period-specific IV validity and the outcome-observation assumption in Equation (S.2) identify the observed average conditional effect for target-population compliers, standardized to $F_{X,hc}^\star$; the required conditional mean equalities then put the observed and counterfactual effects on $\mathcal C_{hck}^\cap$. If no anticipation holds, exclusion holds in the counterfactual state, and $g_{hck}$ is correctly specified,

\[
\beta_{hck}^\star-g_{hck}
=\theta_{hck}^{\cap,\star}.
\tag{S.11}
\]

The complier population in (S.11) can differ across values of $k$, so reading the full event-study path as a path for one fixed complier population requires a stronger construction: the intersection across every reported horizon. Support must hold in the reported period and in every period used to construct $F_{X,hc}^\star$ or $g_{hck}$, and is therefore horizon specific.

Let $\mathcal G_k$ denote the cohorts that contribute at horizon $k$, and $\mathcal H_c$ the event units in cohort $c$. The cohort weights are nonnegative and satisfy $\sum_{c\in\mathcal G_k}\omega_{c\mid k}=1$; within a cohort, $\sum_{h\in\mathcal H_c}\omega_{hc\mid k}=1$. The aggregate is

\[
\theta_k
=\sum_{c\in\mathcal G_k}\omega_{c\mid k}
\sum_{h\in\mathcal H_c}\omega_{hc\mid k}
\left\{
\beta_{hck}^\star
-g_{hck}
\right\}.
\tag{S.12}
\]

There is another possible aggregation, but it answers a different weighted question. If reduced-form and first-stage moments are pooled before the ratio is taken, then for generic components $j$ with prespecified weights $a_j$,

\[
\frac{\sum_j a_j\rho_j}{\sum_j a_j\pi_j}
=\sum_j
\frac{a_j\pi_j}{\sum_\ell a_\ell\pi_\ell}
\frac{\rho_j}{\pi_j}.
\tag{S.13}
\]

Thus the pooled ratio weights components by their first stages rather than using the weights in (S.12). Changes in cohort composition likewise change the target across horizons unless $\mathcal G_k$ and the cohort weights are held fixed across $k$; trimming a component for lack of support changes it too. Report the cohort set, the weights, and the share of the reference distribution retained.

## E. Estimation and inference

### E.1. Finite-strata estimator

The estimator should mirror the estimand: compute a Wald ratio for each period and each stratum with adequate support, then average those ratios using the chosen reference weights.

For finite covariate strata $x\in\mathcal X^\star$, let $\widehat p_x^\star$ be the reference weight and $\widehat g_1$ the estimated counterfactual post-event LATE without the event. This gives the plug-in estimator

\[
\widehat\beta_t^\star
=\sum_{x\in\mathcal X^\star}
\widehat p_x^\star
\frac{\widehat\rho_t(x)}{\widehat\pi_t(x)},
\qquad
\widehat\theta=\widehat\beta_1^\star-\widehat g_1.
\tag{S.14}
\]

In the simplest two-period case, $\widehat g_1=\widehat\beta_0^\star$; other assumptions about the counterfactual may use several pre-event estimates. In either case, apply the support rule and weights exactly as specified by the target population.

The default reference is the pre-event complier distribution. Let $\widehat f_0(x)$ denote the empirical pre-event covariate mass over the supported strata; its estimated weights are

\[
\widehat p_x^\star
=\frac{\widehat\pi_0(x)\widehat f_0(x)}
{\sum_v\widehat\pi_0(v)\widehat f_0(v)}.
\tag{S.15}
\]

For the pre-event component, substituting (S.15) cancels the cell-specific first stages:

\[
\widehat\beta_0^\star
=\frac{\sum_x\widehat f_0(x)\widehat\rho_0(x)}
{\sum_v\widehat f_0(v)\widehat\pi_0(v)}.
\]

The same cancellation does not carry over to the post-event component, which still contains the supported-stratum ratios $\widehat\rho_1(x)/\widehat\pi_1(x)$. Strong-identification inference should re-estimate reference weights determined by the analysis data; Section E.3 takes a different case, in which the weights are fixed and external.

### E.2. Strong-identification inference

A bootstrap for the complete procedure resamples the original independent units and reconstructs every component determined by the data—the reference distribution, generated instrument, cohort stacks, counterfactual LATE, and aggregation weights. Otherwise it misses covariance among the pieces of $\widehat\theta$. Prespecified support can remain fixed; when support is selected from the data, inference must account for that selection and the researcher must say which population remains represented.

Each resample returns the final contrast, since estimating the components separately and treating them as independent discards their covariance. The reconstructed estimates provide the joint bootstrap distribution, and the reported procedure should state how that distribution is turned into a confidence interval. For an event-time path, reconstruct every horizon in the same resample: the horizon-specific distributions give pointwise intervals, while their joint distribution gives a simultaneous band. Analytic variance calculations face the same issue. A sandwich variance needs stacked scores for the exact average-of-ratios functional; averaging the reduced forms and first stages before taking the ratio defines another estimand with another influence function.

### E.3. Weak-instrument inference

The construction below is limited to one pre-event period, one post-event period, fixed finite strata, fixed external reference weights, and a counterfactual post-event LATE equal to the pre-event LATE. Within that scope, it is directly implementable.

Suppose there are $K$ positively weighted strata and $m=4K$ reduced-form and first-stage moments. For each moment $q_j$, construct the interval

\[
\mathcal I_j=
\left[
\widehat q_j-z_{1-\alpha/(2m)}
\widehat{\operatorname{se}}(\widehat q_j),
\quad
\widehat q_j+z_{1-\alpha/(2m)}
\widehat{\operatorname{se}}(\widehat q_j)
\right].
\tag{S.16}
\]

Provided the marginal normal approximations are valid under the declared sampling and clustering design, Bonferroni's inequality gives the resulting rectangle joint coverage of at least $1-\alpha$. Because the standard errors are inputs to this construction, they must already reflect the application's dependence structure; the normal critical value in (S.16) does not provide a few-cluster correction.

For each period-stratum cell, project its reduced-form interval through its first-stage interval. When the first-stage interval excludes zero, the ratio interval is the range of the four endpoint quotients; when it includes zero, take the period-stratum ratio set to be the whole real line. Calling the resulting set $\mathcal B_t(x)$, the confidence set for the event contrast is

\[
\mathrm{CS}_{1-\alpha}
=\sum_x p_x^\star
\{\mathcal B_1(x)-\mathcal B_0(x)\}.
\tag{S.17}
\]

This construction is conservative because it ignores covariance that might tighten the joint set, and weak first stages make the set wider, sometimes unbounded. In the worked files supplied with this guide, which use fixed external weights, the point estimate is 1.00 and the conservative 95 percent set is approximately $[-1.37,3.50]$. The width of this set is the warning: a point estimate that looks clear may still carry little information about the event contrast.

First-stage intervals have a separate interpretive use: a negative point estimate is a warning, though not a test of monotonicity, while an interval spanning zero is weak and inconclusive. If the entire interval lies at or below zero, the data contradict the assumed positive first stage at the reported confidence level; the causal interpretation then stops unless the instrument is recoded, or the design is redefined and defended.

The companion implementation reproduces Equations (S.16)--(S.17), including the worked input and reported result, but its scope does not cover estimated reference weights, data-dependent strata or support, flexible covariates, nonlinear extrapolation of the counterfactual LATE, overlapping cohort stacks, or estimated cohort weights. Such cases call for identification-robust inference tailored to the estimator actually used. If no justified procedure is available, report the reduced forms and first stages, identifying conventional intervals as strong-identification approximations. Anderson--Rubin procedures remain useful for suitable component equations, but valid coverage for their dependent difference requires a joint construction (Anderson and Rubin 1949; Andrews 2022).

### E.4. Dependence and stacked data

Variance estimation has to follow assignment, sampling, event timing, and instrument construction. Repeated observations on the same person share a contribution, as do duplicate appearances created by cohort stacking; keep the original cluster identifiers and rebuild the stacks inside each resample. With few independent clusters, use a procedure justified for the sampling structure at hand.

Few-cluster correction and weak-IV robustness solve different problems: bootstrapping a conventional ratio statistic by cluster does not turn it into a weak-IV procedure, so it remains a strong-identification procedure. Conversely, the projection in Section E.3 relies on credible intervals for the reduced-form and first-stage moments and cannot, by itself, correct few-cluster size distortion.

\Needspace{0.22\textheight}

### E.5. Generated instruments

Suppose first that the scoring rule was estimated in an external training sample. Conditional on that sample and the estimated rule, the analysis may treat the score as fixed, yielding inference conditional on the training sample; unconditional inference must also incorporate uncertainty from estimating the scoring rule. Across periods, keep the rule's input definition, scale, sign, and any cutoff unchanged because these features define the encouragement intervention and therefore the complier population.

If data from the main analysis also enter construction of the score, resample the original independent units and rebuild the training sample, score, cutoff, and assignment rule. The leave-out level should match the source of dependence: leaving out only the focal observation can be too narrow when patients share an ambulance company, market, event unit, or training cluster. Leaving out the relevant observation or cluster prevents mechanical contribution to its own instrument, although independence and exclusion still rest on an institutional argument.

The post-event period may bring new assignment sources or input patterns. Under a prespecified rule, the analysis can classify them as unsupported, map them using supported predetermined features, or examine them separately with an instrument re-estimated after the event. By contrast, retraining the primary score on post-event exposure choices can make the instrument itself respond to the event; for each period and cohort, report score support, shares by instrument value, first stages, and the share of the reference distribution lacking support.

\begingroup
\small
\setlength{\parskip}{0.28em}

## References

Anderson, T. W., and Herman Rubin. 1949. “Estimation of the Parameters of a Single Equation in a Complete System of Stochastic Equations.” *Annals of Mathematical Statistics* 20 (1): 46–63. <https://doi.org/10.1214/aoms/1177730090>.

Andrews, Donald W. K. 2022. “Identification-Robust Subvector Inference.” Cowles Foundation Discussion Paper 2105. Originally issued 2017. <https://cowles.yale.edu/node/139273>.

Callaway, Brantly, and Pedro H. C. Sant'Anna. 2021. “Difference-in-Differences with Multiple Time Periods.” *Journal of Econometrics* 225 (2): 200–230. <https://doi.org/10.1016/j.jeconom.2020.12.001>.

Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. “Bootstrap-Based Improvements for Inference with Clustered Errors.” *Review of Economics and Statistics* 90 (3): 414–427. <https://doi.org/10.1162/rest.90.3.414>.

Doyle, Joseph J., Jr., John A. Graves, Jonathan Gruber, and Samuel A. Kleiner. 2015. “Measuring Returns to Hospital Care: Evidence from Ambulance Referral Patterns.” *Journal of Political Economy* 123 (1): 170–214. <https://doi.org/10.1086/677756>.

de Chaisemartin, Clément, and Xavier D'Haultfœuille. 2018. “Fuzzy Differences-in-Differences.” *Review of Economic Studies* 85 (2): 999–1028. <https://doi.org/10.1093/restud/rdx049>.

Imbens, Guido W., and Joshua D. Angrist. 1994. “Identification and Estimation of Local Average Treatment Effects.” *Econometrica* 62 (2): 467–475. <https://doi.org/10.2307/2951620>.

Konetzka, R. Tamara, Fan Yang, and Rachel M. Werner. 2019. “Use of Instrumental Variables for Endogenous Treatment at the Provider Level.” *Health Economics* 28 (5): 710–716. <https://doi.org/10.1002/hec.3861>.

Miyaji, Sho. 2024. “Instrumented Difference-in-Differences with Heterogeneous Treatment Effects.” RIEB Discussion Paper Series 2024-22. Kobe University. <https://www.rieb.kobe-u.ac.jp/academic/ra/dp/English/dp2024-22.html>.

Rambachan, Ashesh, and Jonathan Roth. 2023. “A More Credible Approach to Parallel Trends.” *Review of Economic Studies* 90 (5): 2555–2591. <https://doi.org/10.1093/restud/rdad018>.

Roth, Jonathan. 2022. “Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends.” *American Economic Review: Insights* 4 (3): 305–322. <https://doi.org/10.1257/aeri.20210236>.

\endgroup
