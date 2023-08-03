# Disruption Networks
This is a public repository for the code associated to the paper: Gerassy-Vainberg et al. A personalized network framework reveals predictive axis of anti-TNF response across diseases.
The analyses done in this study were divided in several reports as follows:  

| Notebook | Input data | Figures |
| --- | --- | --- |
| Data pre-processing | - The microarray raw data generated in this study: GSE186963. <br> - The CyTOF raw data generated in this study: FlowRepository FR-FCM-Z4MQ. <br> - Luminex data – in GitHub data directory. <br> -Additional files for the analysis in the data directory | S2 |
| Dynamics_and_baseline_analyses_DNets | - GEO94648 <br> -Additional input files for the analysis in the data directory | Fig 1, 2, 3, 4a-b <br> S1, S3, S4, S5 |
| scRNAseq_analysis_of_the_response_disrupted_pathways | - The scRNA-seq data generated in this study: PRJNA779701. <br> -Additional input files for the analysis in the data directory | Fig 4c-d, <br> S6 |
| Predictive_signature_and_validation | - GSE20690 <br> - GSE33377 <br> - GSE42296 <br> - in-house CD cohort- qPCR results- in GitHub data directory <br> -scRNAseq processed data <br> - CyTOF processed data <br> - Additional input files for the analysis in the data directory | Fig 4e, Fig 5 <br> S7, S8 |
| Disruption_Networks_functions | | |

<br> 
Whether non-responders' transcriptional profile reflects fundamental routes of IFX resistance, is essential for tailoring treatment. To elucidate molecular mechanisms of individual-specific pathways of treatment non-response, we devised a systematic framework we term ‘Disruption Networks’ which generates a new data-type to provide individual-level information of cell-centered changes in cross-feature relations. The generation of the new data-type relies on studying relations between features across a predefined reference population of individuals (i.e., a population level reference network), and then inferring how these relations differ (i.e., are disrupted) at the single sample level. The new data-type can serve as an input to multiple analyses including integration, differential signal detection, patient stratification based on disruption profile, assessment of disruption in functional modules and evaluation of individual’s molecular network behavior under specific perturbation effects or biological conditions

<br> 
![image](https://github.com/ShiranVaniberg/Disruption_Networks_and_IFX_response/assets/51864609/2536d642-e7f7-4050-a189-112d8461f64a)
