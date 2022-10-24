# changes in version 2.1.1

* Enhancements
  - for platform_design(), one_stage_multiarm() and one_design(), now users can choose to control either
    family-wise or pair-wiser type-I error rate.
  - for platform_design(), the function now returns the number of patients saved in the K+M-experimental     arm trial compared to conducting one K-experimental arm and one M-experimental arm trial separately.
  
* Bug fixes
  - for platform_design(), the error used to occur when K is not equal to M is now fixed.
