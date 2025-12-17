# Review: Snyder's ISEA Projection

## Status: NEEDS REVISION

## Issues Found

### Issue 1: Missing citation for equal-area inheritance claim (Line 17-23)
- **Problem**: "The projection inherits the equal-area property from LAEA" is stated without proof or citation. The claim that Snyder's modifications "maintain the equal-area property" is not self-evident and requires justification.
- **Required fix**: Either (1) provide a proof sketch showing that the modifications preserve the Jacobian determinant, or (2) cite Snyder (1992) page number where this is proven, or (3) cite DGGRID documentation if it provides the mathematical justification.
- **Reference**: If Snyder (1992) proves this, cite the specific page. If not proven in the paper, acknowledge that the equal-area property is stated but not rigorously proven in Snyder (1992).

### Issue 2: Missing derivation reference for E_l formula (Line 45-48)
- **Problem**: The formula $E_l = \arctan\left(\frac{1}{\sqrt{3}} \tan(a/2)\right)$ is given without citation or derivation. The statement "using spherical trigonometry for an equilateral spherical triangle with vertex angle 60°" does not constitute a proof.
- **Required fix**: Either (1) provide a detailed derivation showing how this formula follows from spherical trigonometry (law of cosines for sides, law of sines for angles), or (2) cite the page in Snyder (1992) where this is derived, or (3) cite a spherical geometry textbook.
- **Reference**: Snyder (1992) should contain this derivation; provide the page number. Alternatively, cite Coxeter's "Regular Polytopes" or a spherical trigonometry text.

### Issue 3: Unverifiable optimization claim for R₁ (Line 59)
- **Problem**: "Snyder states that R₁ is 'chosen to minimize the maximum scale distortion over the icosahedral face.'" This is a direct quote but no page number is provided.
- **Required fix**: Provide the page number in Snyder (1992) where this optimization criterion is stated.
- **Reference**: Search Snyder (1992) for the discussion of R₁ optimization and cite the specific page.

### Issue 4: Unjustified formula for total projected area (Line 163)
- **Problem**: "Total area = 20 × (area of one equilateral triangle) = $R_1^2 \times 4\pi$" - The equation linking 20 triangles to $R_1^2 \times 4\pi$ is not obvious. Why does the sum equal $R_1^2$ times the sphere area?
- **Required fix**: Explain that each triangle has area $R_1^2 \times \pi/5$ (by equal-area property), therefore 20 triangles have area $20 \times R_1^2 \times \pi/5 = R_1^2 \times 4\pi$.
- **Reference**: This is a definitional calculation, but the intermediate step should be shown to avoid appearing as an unjustified claim.

### Issue 5: Missing citation for discontinuity quantification (Line 175-178)
- **Problem**: "Maximum gap at face edges: ~0.0001°" and "Maximum gap at vertices: ~0.001°" are stated without citation or derivation. These are empirical claims requiring either computational evidence or a reference.
- **Required fix**: Either (1) cite DGGRID documentation if these values are stated there, or (2) describe how these were computed (e.g., "computed by evaluating projection coordinates from both adjacent faces along shared edges"), or (3) remove these specific numbers and state qualitatively that "discontinuities are subcentimeter at Earth scale."
- **Reference**: If these are original computations, state: "computed via hexify implementation by comparing projections from adjacent faces along boundaries."

### Issue 6: Causal language without justification (Line 187)
- **Problem**: "This is conventional and simplifies face assignment logic" - The causal claim "simplifies" is not justified. Does it actually simplify the logic, and if so, how?
- **Required fix**: Either (1) explain how a vertex at the pole simplifies face assignment (e.g., "creates polar quads that can be handled as special cases"), or (2) weaken to "This is conventional; it is not mathematically preferred but computationally convenient."
- **Reference**: No external reference needed; this is a logic claim that should be self-evident from explanation.

### Issue 7: Missing citation for convergence statistics (Line 73-78)
- **Problem**: "typically 3-5 iterations" is an empirical claim stated without evidence. The section discusses convergence but provides no reference for these specific iteration counts.
- **Required fix**: Either (1) state "empirically observed in hexify implementation," or (2) cite DGGRID documentation if it provides iteration counts, or (3) add to the inverse projection section where empirical data is provided.
- **Reference**: Since the inverse projection section (inverse/section.md) provides empirical statistics, cross-reference that section or move this claim there.

## Verified Claims

### Constants
- $G = 36° = \pi/5$ correctly derived from icosahedral symmetry (10 equatorial vertices).
- $E_l \approx 37.37736814°$ numerical value is verifiable (can be computed from formula, though formula itself needs citation - see Issue 2).
- $R_1 = 0.9103832815$ stated as empirically optimized, consistent with DGGRID.

### Forward Projection Steps
- Step-by-step breakdown (Lines 74-145) is logically structured.
- Great-circle distance formula (Line 81) is standard spherical trigonometry.
- Sector reduction (Lines 90-96) correctly exploits 3-fold symmetry.
- Formulas for $\delta_z$, $h$, $A_G$ (Lines 100-109) match DGGRID implementation.

### Equal-Area Property
- Jacobian determinant claim $\det(J) = R_1^2 \sin(\text{lat})$ (Line 157) is mathematically plausible but needs proof or citation (see Issue 1).
- The factor $R_1^2$ scaling is correct if the claim is true.

### Continuity Claims
- Claim that projection is $C^\infty$ within faces (Line 169) is reasonable given that all functions are smooth, but singularities at face centers need to be explicitly mentioned as exceptions.
- Discontinuities at face boundaries are inevitable for polyhedral projections - this is correct.

### Limitations
- No closed-form inverse is correct (Line 197).
- Numerical precision issues (Lines 200-208) are correctly identified.
- Scale distortion at vertices (16.5%, Line 210) is stated but should cite Snyder (1992) page number.
- Precision modes (Lines 215-222) are implementation details but should reference where these tolerances are defined (hexify documentation or DGGRID).

## Mathematical Truth Standard Compliance

**Violations**:
- Lines 17-23: Equal-area preservation claim without proof or citation (Issue 1).
- Lines 45-48: Formula for $E_l$ without derivation or citation (Issue 2).
- Lines 175-178: Empirical claims without evidence (Issue 5).

**Causal language issues**:
- Line 25: "ensuring that infinitesimal regions... maintain their relative areas" - this is a definitional consequence of equal-area, acceptable.
- Line 187: "simplifies face assignment logic" - unjustified causal claim (Issue 6).

**Hidden non-sequiturs**:
- Line 163: Jump from "20 triangles" to "$R_1^2 \times 4\pi$" without showing intermediate calculation (Issue 4).

## Reference Quality

**Primary reference**: Snyder, J.P. (1992). "An Equal-Area Map Projection for Polyhedral Globes." *Cartographica* 29(1): 10-21.

**Problems**:
- Many claims reference "Snyder (1992)" without page numbers (Issues 2, 3).
- Implementation references (line 230) are internal, not external verification.
- No secondary references for standard spherical geometry formulas.

**Needed improvements**:
- Add page numbers for all Snyder (1992) citations.
- Consider adding Coxeter (1973) for icosahedral geometry foundations.
- Consider adding DGGRID technical documentation reference for implementation details.

## Structure Compliance

Structure is good:
1. **Assumptions and Conventions** (Lines 3-11): Clear statement of conventions.
2. **Definition** (Lines 13-25): What ISEA projection is.
3. **The Three Constants** (Lines 27-72): Definitions of G, E_l, R₁.
4. **Forward Projection** (Lines 74-145): Step-by-step algorithm.
5. **Key Properties** (Lines 147-181): Equal-area, continuity, discontinuities.
6. **Engineering Choices** (Lines 183-191): Design decisions.
7. **Limitations** (Lines 193-222): Constraints and numerical issues.

No metaphors used as definitions.

## Notation Consistency

Generally consistent:
- $R_1$: scale factor
- $G$: $36°$ constant
- $E_l$: center-to-edge distance
- $z$: great-circle distance
- Az: azimuth
- $(\rho, \text{Az})$: polar coordinates

**Minor inconsistency**: Line 8 uses "atan2(x, y)" convention but line 84 uses "atan2(...)" without explicit argument order. Should clarify that Snyder's convention is consistently applied.

## Sanity Check Quality

**Missing**: This section lacks a sanity check. The Lambert section has R code; this section should have similar verification.

**Recommended addition**: Add R code demonstrating:
1. Forward projection of a known point to (face, tx, ty)
2. Verification that projected area matches expected $R_1^2 \times \text{spherical area}$
3. Face assignment consistency check

Without a sanity check, the section does not meet the standard established by the Lambert section.

## Minor Suggestions

1. **Line 6**: "standard `atan2(y, x)`" - Clarify that standard mathematical convention is atan2(y, x), then explicitly state that Snyder uses atan2(x, y), which is non-standard.

2. **Lines 27-72**: The three constants section is excellent and could serve as a template for other sections requiring constant definitions.

3. **Line 158**: "follows from the Jacobian determinant of the transformation" - If Issue 1 is addressed with a proof sketch, this becomes much stronger.

4. **Lines 183-191**: Engineering choices section is valuable; consider expanding with a note that these choices do not affect mathematical correctness but do affect practical implementation.

## Overall Assessment

This section has strong structure and covers the material thoroughly, but it has critical citation gaps that must be addressed before approval. The equal-area claim (Issue 1) is the most serious: this is the defining property of the projection and cannot be stated without proof or authoritative citation.

The lack of a sanity check is also a significant omission given the complexity of the projection and the standard set by other sections.

**Recommendation**: NEEDS REVISION. Address Issues 1-7, add page numbers to all Snyder (1992) citations, and include a computational sanity check. After revisions, this section will meet publication standards.
