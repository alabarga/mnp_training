# mnp_training

Collection of scripts used to train and validate the classifier presented in [DNA methylation-based classification of central nervous system tumours](https://www.nature.com/articles/nature26000).

##### classifier training and cross validation

```
preprocessing.R
```
Reads raw data, performs normalization, basic filtering and batch effect adjustment.
Normalized and filtered data is stored in `./results`

```
training.R
```
Trains the Random Forest classifier and stores the final classifier in `./results`

```
cross_validation.R
```
Performs nested cross-validation and stores the results in `./CV`

```
calibration.R
```
Evaluates the results of the cross validation and fits the calibration model that is stored in `./results`.
Compiles a final report in `CVresults.html`.

##### tumor purity estimation

To train a Random Forest regression model that predicts ABSOLUTE purity of TCGA brain tumor samples, see `purity.html` or 
`purity.Rmd`

