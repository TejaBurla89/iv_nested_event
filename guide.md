---
title: "Instrumental Variables Before and After an Event"
subtitle: "Identification with Endogenous Exposure"
author: "Sriteja Burla"
date: "August 2026"
lang: en-US
---


## Summary

Suppose a hospital is acquired. Did the acquisition change the effect of receiving care there? A direct comparison with patients treated elsewhere is difficult because hospital destination is not random. Ambulance dispatch offers a possible instrument: otherwise similar patients may be assigned to companies whose hospital preferences shift which hospital initially receives them (Doyle et al. 2015).

One could estimate the ambulance IV before and after acquisition and subtract the Wald ratios. Each ratio can identify a local effect for patients whose destination responds to the instrument in that period. To attribute the difference to acquisition, the estimates must refer to the same target population and comparable compliers. We also need to specify how the LATE would have evolved without acquisition, with exclusion holding in that counterfactual state as well as in the observed periods.

The benchmark in this guide assumes that the same units comply before the event and under both post-event states, and that the LATE would have remained at its pre-event value without the event. Under the full conditions stated below, subtracting the standardized pre-event Wald ratio from the standardized post-event ratio then identifies the event effect. Sensitivity analysis allows departures from the assumptions about who complies and what the counterfactual LATE would have been.

## 1. The empirical setting

### 1.1. Roles of the event and instrument

I use a hospital acquired at a known date as the running example: the acquisition is the event, and emergency patients are the units whose hospital exposure may change. Let $D=1$ denote initial care at the focal hospital, while $D=0$ denotes a prespecified alternative exposure—for example, initial care at one comparison hospital or at hospitals drawn from a fixed eligible-hospital distribution. Let $Z\in\{0,1\}$ denote a prespecified encouragement derived from ambulance-company assignment or company preferences.

Ambulance assignment shifts hospital destination within a period. Acquisition may change care at the focal hospital and at the alternatives used for comparison. The role of $Z$ is therefore to identify the effect of hospital exposure in each period. Attributing a change in that effect to acquisition requires assumptions that link the periods.

In fuzzy or instrumented DiD, the instrument comes from a policy or group change (de Chaisemartin and D'Haultfœuille 2018; Miyaji 2024). Here, the instrument operates within each period. The event changes the environment in which the exposure effects are compared.

\Needspace{0.62\textheight}

**Figure 1. The event and the instrument play different causal roles**

![](figures/figure1_nested_event_roles.pdf){width=88%}

*Notes:* The IV shifts exposure within a period. The event may change both the exposure process and the effect of exposure. Exclusion rules out an effect of $Z$ on $Y$ through channels outside the prespecified exposure contrast.

The same problem arises when a workplace policy prompts workers to sort across establishments, when families choose among schools that adopt programs, or when a reform changes how courts assign cases.

### 1.2. Target population and intervention definitions

Define the **target population** by a rule that remains fixed across the event and throughout the analyzed periods. In the hospital application, it might include all eligible emergency patients originating in prespecified pickup areas around the hospital, whether they receive care at the focal hospital or are ultimately taken elsewhere.

Determine eligibility before ambulance assignment. Defining the study population by where patients actually receive care would condition on the exposure that the instrument changes. Predetermined market, cohort, or source-population indicators can refine the target population. Fixed effects for the provider that actually treats a patient raise a different problem: they condition on an endogenous response and may absorb the first stage or select patients on unobserved determinants of outcomes. Konetzka, Yang, and Werner (2019) discuss the related problem of using a patient-level instrument to study a provider-level attribute.

The instrument and exposure must mean the same thing in every period. The two values of $Z$ should represent the same encouragement, and $D=1$ and $D=0$ should describe the same exposure contrast. If $D=0$ pools several hospitals, a change in their shares may change what alternative care means. Fixing the distribution of alternative hospitals addresses this problem. When that is infeasible, researchers must justify treating care at these hospitals as equivalent versions of the exposure, or distinguish among them in the exposure definition.

The data may be a panel or repeated cross sections; eligibility and the availability of outcomes may change over time. State how eligible units enter the records and which outcomes are observed. If exposure affects survival, attrition, or missingness, it also changes whose outcomes contribute to the estimate. That selection requires a separate analysis.

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

The IV conditions must hold in the analyzed sample in every period used in the design. Section 2.2 states the additional requirements when outcomes are selectively observed or when the estimand covers all compliers in the target population. Report the instrument propensity, reduced form, first stage, and evidence on instrument strength for each period. A pooled first stage can conceal a weak component.

### 2.2. Standardization and outcome observation

Raw Wald ratios can change even when the conditional LATEs do not, simply because the observed covariate distribution of compliers changes. Choose a reference distribution $F_X^\star$ before examining outcomes, then average the conditional ratios over the part of its domain with common support:

\[
\beta_t^\star=E_{F_X^\star}[\beta_t(X)],
\qquad
\Delta^{IV,\star}=\beta_1^\star-\beta_0^\star.
\tag{2}
\]

The covariate distribution of pre-event compliers is a useful default. Given the observed pre-event distribution $F_{X0}$, the reference measure is proportional to $\pi_0(x)dF_{X0}(x)$. Restrict it to strata that support both period-specific ratios. Researchers can choose another policy-relevant distribution, but should explain that choice and report any trimming.

Standardization gives the pre- and post-event Wald ratios the same observed covariate distribution. To extend the interpretation to all compliers in the target population, their conditional mean exposure effect must equal that among compliers with observed outcomes. Outcomes must have a positive probability of being observed in every target stratum, and exclusion must hold for the target population as well as for the analyzed records. Selective survival or attrition will often require an explicit model or sensitivity analysis.

### 2.3. Complier populations and the post-event counterfactual

Patients with the same observed covariates may still differ in whether they comply with the instrument. Thus, standardization alone does not ensure that $\Delta^{IV,\star}$ compares the same compliers. The benchmark assumes that the same units would comply before the event and after it under both the observed and counterfactual states. A weaker approach allows membership to change if the relevant conditional mean exposure effects remain equal across the complier populations. Sections B.2 and B.3 state these conditions.

Even with comparable complier populations, a second question remains: what would the post-event LATE have been at the same date without the event? Let $L_1^{E,\cap,\star}$ denote the standardized post-event LATE under the event and $L_1^{0,\cap,\star}$ that missing counterfactual; both refer to units who comply in all three relevant states and use the reference covariate distribution. The event estimand is

\[
\theta^{\cap,\star}=L_1^{E,\cap,\star}-L_1^{0,\cap,\star}.
\tag{3}
\]

After the event, the data reveal only the first term. The simplest assumption sets the missing counterfactual equal to the comparable pre-event LATE. With several pre-event estimates, researchers can instead extrapolate the pre-event LATE path linearly or bound departures from it. The restriction concerns how the LATE would have evolved. A parallel-trends assumption about average untreated outcomes does not supply that restriction.

Exclusion must hold in the counterfactual post-event state too: holding exposure fixed, $Z$ has no effect on the outcome that would be observed without the event. Institutional evidence must support this restriction because the state is unobserved.

No anticipation gives the pre-event estimates their counterfactual meaning. Announcements, preparation, or early changes in ambulance behavior may affect periods close to acquisition. If so, use an earlier reference period or redefine the transition window.

\Needspace{0.58\textheight}

**Figure 2. The observed IV path leaves one counterfactual LATE missing**

![](figures/figure2_missing_late_path.pdf){width=88%}

*Notes:* The observed pre-event and post-event Wald ratios identify period-specific local effects under their respective IV conditions. An event interpretation also requires a comparable complier population and a specification for how the post-event LATE would have evolved without the event. The dashed line is schematic: applications may hold the pre-event LATE fixed, extrapolate its earlier path, or bound plausible counterfactual values.

Under the benchmark conditions above, including a common complier population and a counterfactual post-event LATE equal to the pre-event LATE, the event effect is

\[
\boxed{\theta^{\cap,\star}=\beta_1^\star-\beta_0^\star.}
\tag{4}
\]

A different counterfactual path requires replacing the pre-event ratio with the corresponding counterfactual value. In either case, the interpretation also rests on the stated restrictions on outcome observation, complier membership, anticipation, and exclusion.

### 2.4. What the estimates identify

Table 1 separates the claims supported at each stage of the design.

\Needspace{0.31\textheight}

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

When researchers cannot defend the final row, they can still report a standardized change in period-specific LATEs: exposure effects changed for the complier populations reached by the instrument. Decide before estimation which interpretation the design aims to support and what to report if the stronger assumptions prove unconvincing.

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

Here $\Delta^{0,\star}$ is the change in the standardized LATE that would have occurred without the event for the common complier population. $B^{C,\star}$ measures how changes in who complies affect the comparison. $B^{O,\star}$ measures changes in the gap between observed compliers and compliers in the target population. The decomposition presumes that the period-specific IV conditions and definitions of the interventions are valid.

Historical IV estimates can inform $\Delta^{0,\star}$ if they use the same instrument, exposure, target population, and assumptions about complier comparability. Covariate balance and standardization help assess changes in observed composition. First-stage changes reveal changes in the size of the complier group, but equal first stages can hide substantial entry and exit. Sampling rates, attrition patterns, and checks on which outcomes are observed help assess $B^{O,\star}$.

The three gaps can instead be bounded jointly. If their absolute values are bounded by $h_0$, $h_C$, and $h_O$, then

\[
\theta^{\cap,\star}
\in
\left[\Delta^{IV,\star}-H,\ \Delta^{IV,\star}+H\right],
\qquad H=h_0+h_C+h_O.
\tag{6}
\]

Vary the discrepancies jointly because they may occur together. Rambachan and Roth (2023) provide a template for calibrating restrictions on the counterfactual LATE path. Applying that approach here requires inference for the resulting IV estimand.

\Needspace{0.74\textheight}

### 3.2. Simulation design and results

Figure 3 shows why each assumption matters. The simulation follows four repeated cross sections around one event. The IV is valid in every period, and the event causes a change in LATE of one in every scenario. The scenarios separately vary observed covariate composition, unobserved complier types, and the path the LATE would have followed without the event.

\Needspace{0.56\textheight}

**Figure 3. Each assumption addresses a different source of bias**

![](figures/figure3_standardization_and_continuation.png)

*Notes:* Markers are means across 2,000 replications; vertical lines mark the 2.5th and 97.5th percentiles. Each period has 2,500 observations. The dashed line gives the true change in LATE caused by the event. “Standardized” averages over the pre-event complier covariate distribution. “Pre-event LATE held fixed” uses the comparable pre-event LATE as the counterfactual post-event value; “Linear extrapolation” uses three pre-event LATEs. The replication materials document the data-generating process, fixed seed, and outputs.\footnote{Replication materials: \url{https://github.com/TejaBurla89/iv_nested_event}.}

In the benchmark case, the same units comply in all three states, and all three estimators average about one. When only the observed distribution of $X$ changes, the unstandardized estimate rises to 1.75; standardizing to the pre-event complier distribution brings it back to 1.00.

When unobserved complier types turn over, the standardized estimate remains near 1.80 even though the observed covariate distribution does not change. In the last scenario, the LATE would have risen without the event. Holding the pre-event value fixed gives 1.40; linear extrapolation from three pre-event LATEs gives 0.99.

## 4. Multiple periods and staggered events

### 4.1. Event-time local average treatment effects

With more periods, estimate a path of standardized LATEs by event time. Its pre-event portion shows how the local effect was evolving. Report the first stage and common support alongside these estimates to show how they changed over time.

The pre-event estimates can guide a choice of counterfactual, but cannot establish what would have happened after the event. Pretrend tests may have little power. Choosing an extrapolation because a specification passes such a test also affects subsequent inference (Roth 2022). Specify the reference period, transition window, counterfactual assumption, and sensitivity range in advance.

### 4.2. Cohort-specific estimation and aggregation

With staggered events, the cohorts contributing to an estimate may change across event times. Start by estimating effects separately for each cohort and event time, following the approach to group-time effects in Callaway and Sant'Anna (2021). Within each cohort, estimate conditional reduced forms and first stages for every institution or other event unit. Form Wald ratios in covariate strata with adequate support, standardize them, and then average across event units.

Taking ratios before averaging preserves the chosen weights. Pooling the reduced forms and first stages first gives more weight to components with larger first stages and generally changes the estimand. Use prespecified weights to aggregate the cohort-specific effects. Keeping the same cohorts and weights at each horizon makes comparisons across event time easier to interpret.

Each cohort-period estimate still needs adequate first-stage and covariate support, adequate outcome observation, and credible assumptions relating its complier population to the others. If support disappears, report the share of the reference distribution that is no longer represented. Either justify extrapolating to that part of the population or state that the estimand now covers only the supported set.

Stacking can place the same original observation in several cohort comparisons. Those copies carry the same sampling information, so retain the original cluster identifiers and rebuild the stacks within each bootstrap resample. Assess calendar-time shocks when they change the exposure contrast, compliance behavior, which outcomes are observed, or instrument validity.

## 5. Estimation, inference, and reporting

### 5.1. Estimation

Before looking at outcomes, fix the event date and transition window, and define the target population, instrument, and exposure contrast. Specify the covariates, support rule, and reference distribution. The plan should explain how outcomes are observed, why the complier populations can be compared, and how the counterfactual post-event LATE will be constructed. For a staggered design, also fix the weights within and across cohorts.

Equation (2) determines the calculation: estimate conditional reduced forms and first stages, form Wald ratios where both periods have adequate support, and average them using the reference weights. If implementing this calculation through regression software, reproduce it directly on the same sample with the same weights to check that the software estimates the intended quantity.

Construct a generated instrument from predetermined information that the event cannot affect. Fix the training population and dates, construction rule, sign, cutoff, and observations or clusters to be omitted. Specify how to handle inputs outside the prespecified support. If data from the main analysis help estimate the score, rebuild it along with the rest of the estimator during resampling. Leaving out the focal observation or cluster prevents it from contributing mechanically to its own instrument. Independence and exclusion still require an institutional argument.

\Needspace{0.29\textheight}

### 5.2. A numerical example

Table 2 works through the calculation for two covariate groups with reference weights of 0.60 and 0.40.

**Table 2. From means by instrument value to the standardized contrast**

| Period | Group | Reference weight | Reduced form | First stage | Wald ratio |
|---|---:|---:|---:|---:|---:|
| Pre | A | 0.60 | 0.36 | 0.20 | 1.80 |
| Pre | B | 0.40 | 0.24 | 0.16 | 1.50 |
| Post | A | 0.60 | 0.58 | 0.20 | 2.90 |
| Post | B | 0.40 | 0.40 | 0.16 | 2.50 |

The pre-event standardized LATE is $0.60(1.80)+0.40(1.50)=1.68$; after the event, it is $0.60(2.90)+0.40(2.50)=2.78$. The difference is 1.10. Under the first three rows of Table 1, this is the standardized change in period-specific LATEs. It is an event effect if the final row also holds and the counterfactual post-event LATE equals the pre-event LATE.

### 5.3. Inference

The components of the final event estimate may share observations or depend on the same reference weights, cohort stacks, or generated instrument. Inference must account for their joint uncertainty, including estimation of the counterfactual. When identification is strong, one option is a cluster bootstrap that resamples the original independent units and re-estimates the full contrast, including the reference distribution, cohort stacks, and generated scores. Weak first stages make the ratios unstable and require a different approach.

With weak first stages, inference must remain valid for the final event estimand despite weak identification. Section E.3 gives a conservative confidence-set procedure for one pre-event period, one post-event period, fixed finite strata, fixed external weights, and a counterfactual post-event LATE equal to the pre-event LATE. It first constructs simultaneous intervals for the reduced forms and first stages, then finds the event contrasts consistent with those intervals. The set becomes unbounded when the data cannot exclude a zero denominator.

Estimated weights, support selected from the data, and staggered aggregation fall outside this procedure. They require inference for the full estimator being used. Without such a procedure, report the reduced forms and first stages and label conventional intervals as approximations that require strong identification.

Choose clusters to reflect assignment and dependence, which may occur at the ambulance-company, market, or hospital level, or at more than one of these levels. Clustering only by patient may be insufficient. With few clusters, a wild cluster bootstrap may improve inference for reduced forms and first stages when its conditions hold (Cameron, Gelbach, and Miller 2008). It does not resolve weak identification when applied to a conventional ratio statistic.

### 5.4. Reporting

Report enough information for a reader to reconstruct each event estimate: outcome and exposure means by instrument value, reduced forms, first stages, reference weights, and standardized period-specific ratios. Give observation rates and support losses separately, including the share of the reference distribution retained. Identify the estimates used to construct the counterfactual LATE. A staggered design should also report cohort weights and clustering. For a generated instrument, describe the training sample and state whether inference repeats the steps used to construct it.

Use the interpretation supported by Table 1 when describing the result. With the full identification argument, the estimate is a “change in the standardized LATE caused by the event for units who would comply in all three states.” If the evidence supports only weaker claims about complier comparability or the counterfactual LATE, call it a “standardized change in period-specific LATEs.”

## 6. Conclusion

In the hospital example, a change in the ambulance IV estimate can reflect a change in care, a change in the patients whose destination responds to ambulance assignment, or a change in which outcomes are observed. The LATE may also have changed without acquisition. Attributing the standardized IV difference to acquisition requires the full set of conditions in Table 1 and a counterfactual LATE held at its pre-event value. When the evidence supports a weaker interpretation, report the period-specific LATEs and use sensitivity analysis to show how assumptions about the missing counterfactual and the complier populations affect the event estimate.

\vspace{2.2em}

\noindent\rule{0.34\linewidth}{0.4pt}

\vspace{1.0em}

**Suggested citation**

Burla, Sriteja. 2026. *Instrumental Variables Before and After an Event: Identification with Endogenous Exposure*. Version 1.0, August 2026.

\newpage

# Technical Supplement

## A. Setup and within-period IV

Consider one institution, or another unit at which the event occurs. Let $\mathcal U$ denote its target population. Membership in $\mathcal U$ is determined before the instrument and event. The index $a\in\{0,E\}$ records whether the event has occurred by the post-event date: at $t=1$, $a=E$ is observed and $a=0$ is counterfactual. Under no anticipation (Section B.2), the two paths coincide at $t=0$, so the observed states are $a_0=0$ and $a_1=E$.

For individual $i$ in period $t$, let $D_{it}(a,z)\in\{0,1\}$ denote potential exposure under binary instrument value $z$ and $Y_{it}(a,d)$ the potential outcome under exposure $d$. This notation imposes exclusion: holding exposure fixed, the instrument has no further effect on the outcome. Exclusion is part of the IV assumptions in each observed period. Extending it to the counterfactual post-event state requires a separate assumption.

Keep the instrument coding fixed across periods, with $z=1$ denoting the value intended to increase exposure; monotonicity then requires $D_{it}(a,1)\geq D_{it}(a,0)$. Define the complier set and exposure effect by

\[
\mathcal C_t(a)=\{i\in\mathcal U:D_{it}(a,1)>D_{it}(a,0)\},
\qquad
\tau_{it}(a)=Y_{it}(a,1)-Y_{it}(a,0).
\tag{S.1}
\]

Under the usual binary-IV conditions, Equation (1) identifies the conditional mean effect of exposure among that period's observed compliers (Imbens and Angrist 1994). Section 2.1 discusses how to assess these conditions and report estimates for each period.

The indicator $O_{it}$ records whether an eligible unit's outcome is observed. Conditional on predetermined covariates, the benchmark assumes that neither the instrument nor exposure affects outcome observation. If either affects survival, attrition, or missingness, the analysis needs a separate argument to address selection.

## B. Standardization and event identification

### B.1. Standardization and outcome observation

Let $\mathcal X^\star$ contain the covariate values supported in both observed periods, with reference distribution $F_X^\star$ as defined in Section 2.2.

To extend the result to all compliers in the target population, assume that their conditional mean effect equals that of compliers whose outcomes are observed:

\[
E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i=x,O_{it}=1]
=E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i=x].
\tag{S.2}
\]

Under Equation (S.2), the standardized Wald ratio averages the conditional effects of exposure among all compliers in the target population, using $F_X^\star$:

\[
\beta_t^\star
=E_{F_X^\star}\!\left[
E[\tau_{it}(a_t)\mid i\in\mathcal C_t(a_t),X_i]
\right].
\tag{S.3}
\]

Equation (S.2) restricts only the conditional mean. Equating the joint distributions of potential exposures and outcomes for observed and target-population compliers would require a stronger assumption. Both versions require a positive probability of observing outcomes wherever $F_X^\star$ assigns weight.

### B.2. A common complier population and the post-event counterfactual

Let $\mathcal C^\cap$ be the intersection of the complier sets before the event and under the two possible post-event states:

\[
\mathcal C^\cap
=\mathcal C_0(0)\cap\mathcal C_1(0)\cap\mathcal C_1(E).
\tag{S.4}
\]

The main result assumes that all three sets equal $\mathcal C^\cap$ at every value of $X$ with positive reference weight. For these compliers, define

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

for every relevant $z$ and $d$. These equalities allow the observed pre-event quantities in Equations (S.4)--(S.6) to be written under the no-event state. To specify the missing $L_1^{0,\cap,\star}$, researchers also need a rule that uses the comparable pre-event LATE history. The simplest rule assumes $L_1^{0,\cap,\star}=L_0^{0,\cap,\star}$.

Exclusion must hold in the counterfactual post-event state as well. Because this state is unobserved, the justification must come from institutional evidence.

**Identification result.** Suppose period-specific IV validity and Equation (S.2) hold. If the three complier sets are equal, there is no anticipation, exclusion extends to the counterfactual post-event state, and the counterfactual post-event LATE equals the pre-event LATE, then

\[
\beta_1^\star-\beta_0^\star
=L_1^{E,\cap,\star}-L_0^{0,\cap,\star}
=L_1^{E,\cap,\star}-L_1^{0,\cap,\star}
=\theta^{\cap,\star}.
\tag{S.6}
\]

A different assumption about the counterfactual LATE replaces $\beta_0^\star$ with a prespecified value constructed from the pre-event LATE history.

### B.3. Allowing the complier population to change

Identification can allow the three complier sets to contain different units if the relevant mean effects are equal. To state this weaker condition, write, for any group $G$,

\[
M_{t,a}^{G}(x)
=E[\tau_{it}(a)\mid i\in G,X_i=x].
\]

For this identification argument, impose the following conditional mean equalities:

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

The first equality says that pre-event compliers have the same pre-event mean effect as those who would comply at the post-event date without the event. The counterfactual assumption then specifies how the mean effect for this latter population would evolve. Assuming it remains at its pre-event value gives $M_{1,0}^{\mathcal C_1(0)}(x)=M_{0,0}^{\mathcal C_1(0)}(x)$.

The second equality says that the counterfactual post-event mean effect also holds for $\mathcal C^\cap$. The third imposes the corresponding equality for the post-event mean under the event. These restrictions allow the complier populations to differ in membership.

All three restrictions apply wherever $F_X^\star$ assigns weight, and $\mathcal C^\cap$ must have positive probability there. When group membership differs, equal conditional mean effects need a substantive justification.

## C. Decomposition and sensitivity

### C.1. Decomposition

Let $m_t^O$ denote the conditional effects among observed compliers in period $t$, averaged over the reference distribution. Under period-specific IV validity, $m_t^O=\beta_t^\star$. Let $m_t^U$ denote the corresponding average for all period-$t$ compliers in the target population. The three discrepancies in Equation (5) are

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

Adding and subtracting these effects gives Equation (5). A constant counterfactual LATE for the common complier population sets $\Delta^{0,\star}=0$. Equality of the three complier sets sets $B^{C,\star}=0$, and Equation (S.2) sets $B^{O,\star}=0$.

Under Section B.3's conditional mean equalities and its assumption that the no-event mean effect for $\mathcal C_1(0)$ is unchanged between periods, $\Delta^{0,\star}+B^{C,\star}=0$. Neither term need be zero separately.

### C.2. Sensitivity region

Let $\mathcal H$ contain the prespecified plausible triples $(\delta,b_C,b_O)$ for these discrepancies. The corresponding set of event effects is

\[
\Theta(\mathcal H)
=\left\{
\Delta^{IV,\star}-\delta-b_C-b_O:
(\delta,b_C,b_O)\in\mathcal H
\right\}.
\tag{S.9}
\]

Equation (6) uses symmetric rectangular bounds. More generally, $\mathcal H$ can allow asymmetric bounds, impose signs, or restrict how the discrepancies vary together. Historical IV estimates can inform the counterfactual LATE component when they use the same instrument, exposure, target population, reference distribution, and assumptions about complier comparability.

## D. Multiple periods and staggered events

With multiple periods, apply the two-period argument separately within each cohort at each event time, then aggregate the estimates.

Let $c$ index event cohorts, $h$ institutions or other event units, and $k$ event time. Write the conditional reduced form and first stage for each component as $\rho_{hck}(x)$ and $\pi_{hck}(x)$. Using reference distribution $F_{X,hc}^\star$, its standardized LATE in the observed period is

\[
\beta_{hck}^\star
=E_{F_{X,hc}^\star}\!\left[
\frac{\rho_{hck}(X)}{\pi_{hck}(X)}
\right].
\tag{S.10}
\]

For a reported post-event horizon, define $\mathcal C_{hck}^\cap=\mathcal C_{hck}(0)\cap\mathcal C_{hck}(E)$. Construct the prespecified counterfactual LATE, $g_{hck}$, from comparable pre-event LATEs for the same complier population and reference distribution. Period-specific IV validity and Equation (S.2) identify the conditional effects among target-population compliers in the observed state, averaged over $F_{X,hc}^\star$. The required conditional mean equalities allow both the observed and counterfactual effects to apply to $\mathcal C_{hck}^\cap$. If there is no anticipation, exclusion holds in the counterfactual state, and $g_{hck}$ is correctly specified,

\[
\beta_{hck}^\star-g_{hck}
=\theta_{hck}^{\cap,\star}.
\tag{S.11}
\]

The complier population in (S.11) can differ across values of $k$. To interpret the full event-study path for a fixed complier population, use the intersection across every reported horizon. Support must hold in the reported period and in each period used to construct $F_{X,hc}^\star$ or $g_{hck}$. It therefore has to be assessed separately for each horizon.

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

Pooling reduced forms and first stages before taking their ratio changes the weights. For generic components $j$ with prespecified weights $a_j$,

\[
\frac{\sum_j a_j\rho_j}{\sum_j a_j\pi_j}
=\sum_j
\frac{a_j\pi_j}{\sum_\ell a_\ell\pi_\ell}
\frac{\rho_j}{\pi_j}.
\tag{S.13}
\]

The pooled ratio therefore weights each component according to its first stage, changing the weights in (S.12). The target also changes across horizons unless $\mathcal G_k$ and the cohort weights remain fixed across $k$. Dropping a component that lacks support changes the target as well. Report which cohorts contribute, their weights, and the share of the reference distribution retained.

## E. Estimation and inference

### E.1. Finite-strata estimator

Compute a Wald ratio for each period in each stratum with adequate support, then average the ratios using the chosen reference weights.

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

In the simplest two-period case, $\widehat g_1=\widehat\beta_0^\star$. Other counterfactual assumptions may use several pre-event estimates. Apply the support rule and weights specified for the target population.

\Needspace{0.14\textheight}

The default reference is the pre-event complier distribution. Let $\widehat f_0(x)$ denote the empirical pre-event covariate mass over the supported strata. The estimated reference weights are

\[
\widehat p_x^\star
=\frac{\widehat\pi_0(x)\widehat f_0(x)}
{\sum_v\widehat\pi_0(v)\widehat f_0(v)}.
\tag{S.15}
\]

For the pre-event component, substituting (S.15) cancels the first stages within each stratum:

\[
\widehat\beta_0^\star
=\frac{\sum_x\widehat f_0(x)\widehat\rho_0(x)}
{\sum_v\widehat f_0(v)\widehat\pi_0(v)}.
\]

The post-event component still contains the ratios $\widehat\rho_1(x)/\widehat\pi_1(x)$ in the supported strata; its first stages do not cancel. When the reference weights come from the analysis data, inference under strong identification should re-estimate them. Section E.3 considers fixed weights obtained from an external source.

### E.2. Strong-identification inference

Resample the original independent units and reconstruct every component determined by the data: the reference distribution, generated instrument, cohort stacks, counterfactual LATE, and aggregation weights. This preserves covariance among the components of $\widehat\theta$. Prespecified support can remain fixed. If support is selected from the data, inference must account for that selection, and researchers must state which population the estimates represent.

Calculate the final contrast in each resample and state how its bootstrap distribution is used to construct a confidence interval. For an event-time path, reconstruct all horizons within each resample. Use the distribution at each horizon for pointwise intervals and their joint distribution for a simultaneous band. Analytic variance calculations must also retain covariance among the components. A sandwich variance requires stacked scores for the estimator that averages the conditional ratios. Pooling the reduced forms and first stages before taking their ratio changes both the estimand and its influence function.

### E.3. Weak-instrument inference

The following construction can be implemented in the limited case of one pre-event period, one post-event period, fixed finite strata, fixed external reference weights, and a counterfactual post-event LATE equal to the pre-event LATE.

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

If the marginal normal approximations are valid under the stated sampling and clustering design, Bonferroni's inequality gives this rectangle joint coverage of at least $1-\alpha$. The standard errors must account for the application's dependence structure. The normal critical value in (S.16) does not correct inference with few clusters.

\Needspace{0.21\textheight}

Within each period and stratum, form a set of ratios from the reduced-form and first-stage intervals. If the first-stage interval excludes zero, calculate the four quotients of their endpoints; the smallest and largest quotients bound the ratio interval. If the first-stage interval includes zero, use the whole real line. Denote the resulting set by $\mathcal B_t(x)$. The confidence set for the event contrast is

\[
\mathrm{CS}_{1-\alpha}
=\sum_x p_x^\star
\{\mathcal B_1(x)-\mathcal B_0(x)\}.
\tag{S.17}
\]

This construction is conservative because it ignores covariance that could narrow the joint set. Weak first stages widen the set and can make it unbounded. The worked files supplied with this guide use fixed external weights and give a point estimate of 1.00, with a conservative 95 percent set of approximately $[-1.37,3.50]$. Despite the point estimate, this wide set leaves substantial uncertainty about the event contrast.

First-stage intervals also help assess the design. A negative point estimate warrants attention, although it is not a test of monotonicity. An interval spanning zero provides weak and inconclusive evidence. If the entire interval lies at or below zero, the data contradict the assumed positive first stage at the reported confidence level. A causal interpretation then requires recoding the instrument or redefining and defending the design.

The companion implementation reproduces Equations (S.16)--(S.17), including the worked input and reported result. It does not cover estimated reference weights, strata or support chosen from the data, flexible covariates, nonlinear extrapolation of the counterfactual LATE, overlapping cohort stacks, or estimated cohort weights. These cases require identification-robust inference for the estimator used. If no justified procedure is available, report the reduced forms and first stages and label conventional intervals as strong-identification approximations. Anderson--Rubin procedures remain useful for suitable component equations. To obtain valid coverage for a difference between dependent components, however, the procedure must handle them jointly (Anderson and Rubin 1949; Andrews 2022).

### E.4. Dependence and stacked data

Variance estimation must account for assignment, sampling, event timing, and instrument construction. Repeated observations on the same person and copies of an observation in different cohort stacks do not provide independent information. Keep the original cluster identifiers and rebuild the stacks in each resample. With few independent clusters, use a procedure justified for the sampling structure.

A cluster bootstrap of a conventional ratio statistic still requires strong identification. The projection in Section E.3 allows weak instruments but relies on credible intervals for the reduced-form and first-stage moments. It cannot correct size distortions caused by few clusters. The two problems require separate treatment.

\Needspace{0.22\textheight}

### E.5. Generated instruments

Suppose the scoring rule was estimated in an external training sample. The analysis can condition on that sample and the estimated rule and treat the score as fixed. Unconditional inference must also account for uncertainty in the estimated scoring rule. Keep the rule's input definition, scale, sign, and any cutoff unchanged across periods: these features define the encouragement intervention and therefore the complier population.

If the score also uses data from the main analysis, resample the original independent units and rebuild the training sample, score, cutoff, and assignment rule. Choose which observations to leave out according to the source of dependence. Omitting only the focal observation may be insufficient when patients share an ambulance company, market, event unit, or training cluster. Leaving out the relevant observation or cluster prevents it from mechanically contributing to its own instrument. Independence and exclusion still require an institutional argument.

\Needspace{0.16\textheight}

New assignment sources or input patterns may appear after the event. A prespecified rule can classify them as unsupported, map them using supported predetermined features, or examine them separately with an instrument re-estimated after the event. Retraining the primary score on post-event exposure choices can make the instrument itself respond to the event. For each period and cohort, report score support, shares by instrument value, first stages, and the share of the reference distribution lacking support.

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
