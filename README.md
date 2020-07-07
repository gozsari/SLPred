# SLPred: a multi-view subcellular localization prediction tool for human proteins
* SLPred is a multi-view subcellular localization prediction tool for human proteins.
* The tool consists of nine independently developed model for the proteins which have annotation with nine subcellular locations as shown in the subcellular locations figure.
* SLPred exploits the features of nineteen different protein descriptors from the publicly available tools: POSSUM, iFeature and SPMAP.
* Support Vector Machine (SVM) is used to construct probabilistic prediction models, which produces probabilistic scores indicating the localization probability for a query protein sequence. 
* A weighted score is calculated based on the obtained probabilistic scores from seven feature-based probabilistic prediction models (SVMs) by employing weighted mean voting.
* Binary prediction is given by applying thresholding on the weighted score.

* The following figure shows the proposed method
![alt text](https://github.com/gozsari/SLPred/blob/master/images/model_architecture.png)

* The following figure shows the subcellular location and their groups.
* This mapping is formed by considering is a and part of relations in the subcellular location hierarchy.
![alt text](https://github.com/gozsari/SLPred/blob/master/images/subcellular_locations.png)
