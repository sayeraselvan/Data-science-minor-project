# Modelling the impact of DNA sequence repeats on gene inversions and gene-strand bias in bacteria


[![Version](https://img.shields.io/badge/version-1.0.0-blue)](https://github.com/sayeraselvan/Data-science-minor-project)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/sayeraselvan/Data-science-minor-project/blob/main/LICENSE.md)

A computational model simulating how DNA sequence repeats influence gene inversions and strand bias in bacterial genomes.

## Abstract
This project investigates how different DNA repeat configurations (direct, inter-inverted, intra-inverted) affect gene-strand bias through inversion events. The model demonstrates:
- **Null Model Behavior**: Without selection pressure, gene-strand bias stabilizes at 50-50 distribution
- **Selection Impact**: Imposing constraints through inversion disparity scores maintains existing strand bias
- **Evolutionary Insights**: Provides framework for understanding bacterial adaptation through chromosomal rearrangements 

## Key Features
- Circular bacterial genome simulation with configurable parameters:
  - Ori/Ter positions
  - Repeat distributions (direct, inter/intra-inverted)
  - Gene allocation between strands
- Two operational modes:
  - **Null Model**: Unconstrained inversions
  - **Selection Model**:
    - Gene count imbalance thresholds
    - Inversion Disparity Score (IDS) calculations
    - Penalty limit constraints
- Analytical tools for:
  - Strand bias quantification
  - Inversion event tracking
  - Fitness impact assessment



**Key Configurable Parameters** (via `config.yaml`):
- `inversion_disparity_limit`: 10/25/no-limit
- `repeat_distribution`: Random/Clustered
- `fitness_function`: Normal/Exponential distribution

## Results Interpretation
### Null Model Dynamics 
- Rapid convergence to 50-50 strand distribution
- Demonstrates baseline evolutionary pressure from unconstrained inversions

### Selection Model Outcomes (Figs 5-7)
| Penalty Limit | Strand Bias Stability | Inversion Frequency |
|---------------|-----------------------|---------------------|
| 10            | High                  | Low                 |
| 25            | Moderate              | Medium              |
| no-limit      | Low                   | High                |

**Critical Finding**: Higher penalty limits permit larger inversions while maintaining ancestral strand bias patterns 

## Future Directions
Planned model extensions:
1. GC skewness integration
2. Transcription-replication collision modeling
3. Horizontal gene transfer simulations


## License
MIT License - see [LICENSE](https://github.com/sayeraselvan/Data-science-minor-project/blob/main/LICENSE.md) for details.

