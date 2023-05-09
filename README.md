# SABER-PLISH_barcodes

This repository contains all the scripts used in the computational design of DNA barcode libraries for use in SABER-PLISH.

Documentation to improve its reusability has yet to be completed.

The DNA barcode libraries generated are not uploaded here due to large file sizes. They will instead be uploaded to a private Amazon S3 bucket.

# Using workflow.sh to generate (1) nanobody, (2) probe and (3) implied barcode libraries

1) In the commandline, run `bash workflow.sh`.
2) A prompt requesting length of sequence (L) will be shown. Note that L must be an even number.
3) A prompt requesting degree of orthogonality (k) will be shown. It is recommended that k <= 0.4L - 0.5L.
4) A prompt requesting a boolean to set RCfree argument to True or False will be shown. For high values of L or k, you are recommended to set this to False (e.g. L=42 and k=13 with RCfree=True is too computationally expensive to compute).
5) 
