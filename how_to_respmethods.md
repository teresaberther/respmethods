This is a step-by-step, more conceptual introduction into the code for the robust pipeline for respiration phase-locked analyses of neural and behavioral data provided in this repository. Note this is an exemplary workflow, the code can be adapted for any event-based or continuous analysis of choice, as long as behavioral or neural data was recorded simultaneously with respiration. 
## Event-based analysis pipeline
We start this pipeline based on the following input: 
- continuous, raw respiration time series vector
- events structure that contains the timestamp of interest in each trial in samples from start of recording and the corresponding outcome value of interest for each trial 
This exemplary pipeline uses hit rate in a visual near-threshold task fluctuating across the respiratory cycle as an example for a behavioral outcome variable of interest. <br>

**Relevant Scripts**: ```Tutorial_DataPrep``` (MATLAB) / ```Tutorial_DataPrepMultiprocessing``` (Python) and ```Tutorial_DataClust``` (both)

### Step 1: Extract respiratory phase 
Starting with the ```Tutorial_DataPrep```script, we load, downsample, smooth and z-score the respiratory time series. The original sampling frequency as well as the sampling frequency after downsampling are set as parameters in the ```config```section in the script head. Then, we extract continuous respiratory phase using the two-point-interpolation function ```two_point_interp``` in the ```Tutorial_DataPrep``` script. The function requires no input besides the z-scored respiratory time series and can be flexibly interchanged with any of the other provided phase extraction methods using the same syntax. 
**Exemplary function call** (equivalent for both MATLAB and python): ```respphase = two_point_interp_(zresp)``` with ```zresp```being the z-scored, smoothed original respiration trace vector and the resulting ```respphase```vector being the extracted phase

**Sanity check:** Plot the z-scored time series overlaid with the phase time series resulting from the phase extraction function.

**Additional tip**: All interpolation functions are dependent on the criterion for peak detection - try changing the criterion value and plotting the z-scored time series with the detected peaks for a sanity check if phase extraction fails with the default parameters. 

### Step 2: Generate surrogate phase vectors using IAAFT 
Next in ```Tutorial_DataPrep```, a number of surrogate time series are computed using the ```generate_surrogates_iaaft```function. The number of surrogates can be set via the ```niter```parameter in the script ```config```head. We recommend a minimum number of 1000 surrogates and use a default setting of 5000 surrogate iterations in the example pipeline. ```generate_surrogates_iaaft```takes the amplitude and frequency distribution of the z-scored respiration time series and generates a phase-shuffled 'mock' signal keeping all other characteristics of the original time series (including the autocorrelation) intact. The function only requires the **z-scored original time series vector** as an input. As our aim with surrogate generation is to create a distribution of signals whose phase is independent from the phase of the original signal, we generate one iaaft solution at a time, extract its phase with the same method used on the empirical phase extraction, and check the surrogate phase congruence with the original time series. The congruence constraint criterion ```plvcrit```can be adjusted in the ```config``` head. 

**Exemplary function calls**:

MATLAB
```matlab
surr = generate_surrogate_iaaft(zresp);
```
and python
```python
surr = iaaft(zresp)
```
with the resulting phase-shuffled surrogate time series ```surr```.

### Step 3: Moving window phase binning 
Last in ```Tutorial_DataPrep```, we use a moving window approach to average the outcome resolved by respiration phase. Using the extracted phase vectors, we determine the concurrent respiratory phase in each trial with the event of interest happened (here: stimulus presentation). Then, we iterate across phases, finding trials whose phase falls into the current iteration phase bin and calculating the mean across all outcome values of these trials. The granularity and overlap of the analysis are determined by the number of bins ```nbin```and bin width ```phw```, respectively, which can be set in the ```config```head of the script. Refer to the accompanying tutorial paper for a detailed description on how to set the binning parameters. The binning is run on both the empirical phase vector as well as all surrogate iterations. 

**Additional tip**: If outliers in the outcome variable across trials are a concern, switch out the aggregation of values for each phase bin using means for a more robust measure, e.g. trimmed means. 

### Step 4: Circular Cluster Permutation testing 
After the phase binning of the outcome variable, we test the significance of phase-modulated effects using our novel cluster permutation approach for circular data as illustrated in ```Tutorial_DataClust```. T-values are computed on the combined empirical and surrogate distributions and piped into the permutation testing function ```CircPerm```. ```CircPerm```requires the empirical and surrogate t-values as separate input arrays, as well as the vector containing the phase bin centres chosen in Step #3. Additionally, the tail(s) and direction of the test can be set with the last parameter, allowing for two sided (pass ```"two_sided"```) or one sided (pass ```"greater"```or ```"lesser"```) testing. 
The function returns a structure ```summary```containing the clustering results as well as the boundaries ```bounds```used for defining significant data points in each phase bin. 
In ```summary```, each row represents a separate cluster and contains four fields: 1. ```.idxs```with the phase bin indices for each cluster, 2. ```.clustmass_stat```with the clustering statistic for each cluster (the sum of z-scores for all points in the cluster), 3. ```p_ecdf```with the empirical cumulative distribution function (ecdf) value for each cluster statistc, and 4. ```.p``` the permutation p-value.  

**Exemplary function calls**:

MATLAB
```matlab
[summary, bounds] = CircPerm(emp_t, perm_t, phasevec, "two_sided");
```
and python
```python
[summary, bounds] = circ_perm(emp_t, perm_t, phasevec, "two_sided")
```

**Additional tip #1**: Note that the smallest attainable p-value in the permutation test is determined by the number of surrogates tested and becomes smaller as the number of surrogates increases - consider this when setting the ```niter``` parameter in Step #2.

**Additional tip #2**: We run the circular cluster permutation test on t-values because they are widely used, easily interpretable, and account for sample size. However, they are based on sample means and SD and can be less robust if outliers or skewness are present in a sample. Thus, we have added some alternate code using a more robust scaling on sample level based on the interquartile range (IQR), as our permutation logic works on values subjected to any scaling as long as the scaling is applied to the combined empirical and surrogate distribution. Note that interpretation differs for classical t-values and robustly scaled values derived from the IQR scaling. 


## Continuous analysis pipeline 
We start this pipeline based on the following input:
- continuous, raw respiration time series vector
- continuous, clean sensor or source-level M/EEG data 

This exemplary pipeline uses global field power in neural resting state data fluctuating across the respiratory cycle as an example for a neural outcome variable of interest, computing phase-amplitude coupling.
Note that preprocessing/cleaning of continuous neural data is not the focus of this tutorial, thus we kindly refer the reader to the relevant [*fieldtrip*](https://www.fieldtriptoolbox.org/tutorial/preproc/continuous/) and [*MNE*](https://mne.tools/stable/auto_tutorials/preprocessing/index.html) tutorials for this step.

**Relevant Scripts**: ```Tutorial_MI```

### Step 1: Extract respiratory phase & Step 2: Generate surrogate phase vectors using IAAFT 
The ```Tutorial_MI```script starts with the exact same two steps as the event-based analysis: continuous phase extraction using interpolation and surrogate time series generation using IAAFT. For details on function calls and parameters, refer to the descriptions provided in the event-based analysis section above. 

### Step 3: Modulation Index calculation 
After extracting phase and generating surrogate time series, we calculate the Modulation Index as introduced by Tort (2010). In our example, we analyse the phase-amplitude coupling of respiration to global field power between 0.5 and 40Hz, but the same logic can also be applied to the amplitude envelope of each channel/sensor/parcel/voxel of sensor/channel-level or source-level data. The MI is computed for both empirical and surrogate respiration phase time series, and the resulting empirical MI values are expressed in units of standard deviations based on the surrogate MI distribution. Similar to the parameter for the event-based analysis, the  ```config```head of the ```Tutorial_MI```script also sets a ```nbin```number of bins to determine the granularity for the MI calculation. We propose a default of 60 bins, but refer the reader to the [original MI publication](https://journals.physiology.org/doi/full/10.1152/jn.00106.2010?iss=4) for more detailed recommendations and considerations on ideal granularity settings.


