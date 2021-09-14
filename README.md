# A Geometric Perspective On Functional Novelty Detection 

Code and data to reproduce the results presented in the paper *A Geometric Perspective On Functional Novelty Detection* by Herrmann and Scheipl.
The full paper can be reproduced with the 00_paper.Rmd in the vignettes folder.

Some of the examples and results are computed directly in the Rmd files in folder vignettes, others are (partly) outsourced to specific R files in folder R:
  - examples in section 2
    - 03_theory.Rmd
  - qualitative experiments in section 3.1
    - 04_qual_analysis.Rmd 
    - exp_real_qual.R
  - quantitative experiments in section 3.2
    - 05_quant_experiments.Rmd
    - exp_quant_complex_batch.R (results depicted in Figure 9) 
    - exp_l2vsl10_batch.R (results depicted in Figure 10 A)
    - exp_arrowhead_batch.R (results depicted in Figure 10 B)
    - exp_quant_complex_ui_batch.R (results depicted in Figure 11)
  - appendix
    - 07_appendix.Rmd
    - exp_sensitivity_extended.R (Appendix B)
    - exp_fdaoutlier_batch.R (Appendix C)
    - exp_outliergram_msplot.R (Appendix D)

File names ending on *batch.R* indicate that the experiments where conducted using package `batchtools` and take some time to compute. 




