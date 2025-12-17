# Theory Section Rewrite - Orchestration Plan

## Objective
Completely rewrite the mathematical foundations vignette with:
- Mathematically valid statements (every claim is definition, theorem with proof sketch, or labeled engineering choice)
- Proper citations for all claims
- Clean, consistent illustrations following strict quality gates
- No preserved explanations unless verified against references

## Agents

| Agent | Topic | Status |
|-------|-------|--------|
| A | Lambert Azimuthal Equal-Area Projection | Research phase |
| B | Snyder ISEA Polyhedral Projection | Research phase |
| C | Icosahedron Geometry and Face Assignment | Research phase |
| D | Apertures and Subdivision Patterns | Research phase |
| E | Cell Indexing and Hierarchical IDs | Research phase |
| F | Inverse Projection and Newton-Raphson | Research phase |
| G | Editor and Sanity Checker | Pending (starts after Phase 1) |

## Phases

### Phase 1: Research (Current)
Each agent gathers references and creates documentation in `01_references/`.

### Phase 2: Draft Writing
Each agent writes their section in `02_drafts/` as standalone markdown.

### Phase 3: Figure Creation
Agents produce code and exports in `03_figures/` with quality gate specs.

### Phase 4: Internal Review
Agent G reviews all drafts and figures, writes review notes in `04_review/`.
Iterate until all sections pass.

### Phase 5: Integration
Merge accepted drafts into `05_final/theory.md`.
Create `05_final/ERRATA.md` documenting what was wrong in original.

## Quality Gates

### Mathematical Claims
- [ ] Every claim is: definition, theorem with proof sketch, or cited
- [ ] No "therefore" without logical implication from stated assumptions
- [ ] No "it can be shown" without proof sketch or citation

### Figures
- [ ] Reproducible from code
- [ ] Exported as SVG + PNG
- [ ] No unintended overlaps
- [ ] Minimal palette: grayscale + 1 accent (+ 1 more only if encoding different concept)
- [ ] Hexagon figures: no self-intersection, consistent spacing/orientation

### References
- [ ] Primary sources preferred
- [ ] Page numbers for books/papers
- [ ] Stable links (not Wikipedia as primary)

## Timeline
Started: 2025-12-17
