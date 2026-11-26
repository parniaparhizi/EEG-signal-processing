
# EEG-Based Analysis of Olfactory Responses in Alzheimer's Disease

### Time-Frequency Analysis of EEG Signals Using MATLAB and EEGLAB

## Overview

This project investigates EEG signal responses to olfactory stimulation across different cognitive conditions, with a focus on potential electrophysiological markers associated with Alzheimer's disease.

The study applies signal processing techniques to examine changes in neural oscillatory activity during exposure to different odors.

The primary focus of this phase is time-frequency power analysis using the Short-Time Fourier Transform (STFT), particularly in the theta and gamma frequency bands.

The analysis provides a foundation for future investigation of phase-amplitude coupling (PAC) as a potential biomarker of cognitive decline.

**Course:** Signals and Systems  
**Institution:** Sharif University of Technology  
**Instructor:** Prof. Hamid Aghajan  
**Semester:** Spring 2025

---

## Research Motivation

Alzheimer's disease is a progressive neurological disorder that affects memory and cognitive function.

Olfactory dysfunction may appear during early stages of cognitive decline because brain regions involved in odor processing are closely connected to memory-related neural networks.

EEG provides a non-invasive method for investigating changes in neural activity during sensory stimulation.

This project explores whether time-varying EEG power patterns differ between healthy individuals and subjects experiencing cognitive impairment.

## Project Objectives

The main objectives include:

- Preprocessing EEG recordings using MATLAB and EEGLAB.
- Removing noise and artifacts from EEG signals.
- Extracting event-related EEG responses to olfactory stimuli.
- Applying time-frequency analysis using STFT.
- Investigating theta and gamma oscillatory activity.
- Comparing neural responses to chocolate and rose odors.
- Examining differences between healthy control, MCI, and Alzheimer's subjects.
- Establishing a foundation for future phase-amplitude coupling analysis.

## Dataset

The analysis uses EEG recordings collected during an olfactory stimulation experiment.

### Subject Groups

Three recordings representing different cognitive conditions were analyzed:

| Group | Description |
|-------|-------------|
| HC | Healthy Control |
| MCI | Mild Cognitive Impairment |
| AD | Mild Alzheimer's Disease |

One subject from each group was included in this phase.

### Experimental Setup

| Parameter | Value |
|-----------|-------|
| EEG Channels | 19 |
| Event Tagging Channel | 1 |
| Sampling Frequency | 250 Hz |
| Odor Conditions | Chocolate and Rose |
| Odor Presentation | 5 seconds |
| Rest Interval | 6 seconds |
| Trials per Odor | 22 |
| Chocolate Event Tag | 5 |
| Rose Event Tag | 6 |

The EEG recordings were synchronized with odor presentation events to enable stimulus-related analysis.

## EEG Preprocessing

EEG preprocessing was performed using MATLAB and EEGLAB.

The preprocessing pipeline included the following steps:

### 1. Data Import

Raw EEG recordings were imported from MATLAB data files.

Standard electrode locations were assigned to the EEG channels.

### 2. Event Extraction

Odor presentation events were identified using the dedicated tagging channel.

Event markers were used to align the EEG recordings with stimulus onset.

### 3. Bandpass Filtering

A bandpass filter was applied to retain the relevant EEG frequency components.

**Frequency range: 1–90 Hz**

### 4. Notch Filtering

A notch filter was applied to suppress power-line interference.

**Frequency range: 48–52 Hz**

### 5. Noisy Channel Detection

Noisy and flatline channels were identified using EEGLAB preprocessing tools.

### 6. Channel Interpolation

Removed channels were reconstructed using spatial interpolation to preserve the electrode montage.

### 7. Average Re-Referencing

EEG recordings were re-referenced using the average reference method.

### 8. Artifact Removal

Artifact Subspace Reconstruction (ASR) and Independent Component Analysis (ICA) were used to identify and reduce unwanted components.

ICLabel was used to assist in component classification.

### 9. Epoch Extraction

EEG signals were segmented around odor presentation events.

The analysis considered time intervals surrounding stimulus onset to investigate event-related neural activity.

### 10. Trial Rejection

Noisy trials were identified and rejected to improve the reliability of subsequent analysis.

---

## Time-Frequency Analysis

### Short-Time Fourier Transform

The Short-Time Fourier Transform (STFT) was used to investigate how EEG frequency content changes over time.

Unlike the conventional Fourier Transform, STFT preserves information about the temporal evolution of different frequency components.

The STFT is defined as:

$$
X(m,\omega)=\sum_n x[n]w[n-m]e^{-j\omega n}
$$

where:

- `x[n]` represents the EEG signal.
- `w[n-m]` represents the analysis window.
- `m` represents the window position.
- `ω` represents angular frequency.

The time-frequency power was calculated as:

$$
P(t,f)=|STFT(t,f)|^2
$$

This representation was used to investigate changes in neural activity during olfactory stimulation.

### Frequency Bands

The analysis focused primarily on two EEG frequency bands:

| Frequency Band | Range |
|----------------|-------|
| Theta | 4–8 Hz |
| Gamma | 30–50 Hz |

Theta oscillations are associated with cognitive processes such as attention and memory.

Gamma oscillations are involved in sensory processing and neural information integration.

### STFT Parameters

The analysis considered:

- Window lengths between 0.5 and 1 second.
- High window overlap (85–95%).
- Time-resolved power estimation.
- Trial-based power analysis.

These parameters provide a trade-off between temporal and frequency resolution.

---

## Power Analysis

Time-varying power was extracted from the STFT representation for the selected frequency bands.

For each condition, power values were organized across time and trials.

Average power was then calculated across trials to investigate the characteristic response of each subject.

The analysis considered:

- Subject group
- Odor type
- Frequency band
- Time relative to stimulus onset

## Experimental Findings

### 1. Healthy Control vs. Cognitive Impairment

The healthy control subject generally exhibited stronger and earlier EEG power responses.

The MCI and Alzheimer's subjects demonstrated weaker or delayed responses, particularly in the theta band.

The Alzheimer's subject showed the most pronounced reduction in oscillatory power.

### 2. Chocolate vs. Rose Stimulation

Chocolate stimulation generally produced stronger power responses than rose stimulation in the analyzed recordings.

This difference was particularly noticeable in the theta and gamma bands.

### 3. Gamma-Band Activity

Gamma-band responses appeared more strongly reduced in the Alzheimer's subject.

This observation suggests that high-frequency neural activity may be affected by cognitive impairment.

### 4. Temporal Response Characteristics

The analysis revealed differences in both response magnitude and timing.

These characteristics may provide useful information for future studies of stimulus-related neural dynamics.

**Important:** These observations are exploratory and are based on one subject per group. They do not establish population-level differences or a validated diagnostic biomarker.

---

## Connection to Phase-Amplitude Coupling

Phase-amplitude coupling (PAC) describes interactions between the phase of low-frequency oscillations and the amplitude of high-frequency oscillations.

For example, theta-gamma coupling investigates whether gamma-band amplitude is modulated by the phase of theta oscillations.

PAC is relevant to neural communication and memory-related processing.

The current project establishes a foundation for PAC analysis by examining the time-varying power characteristics of theta and gamma oscillations.

**Direct PAC estimation is outside the scope of this phase.**

## Tools and Technologies

- MATLAB
- EEGLAB
- Digital Signal Processing
- Short-Time Fourier Transform (STFT)
- EEG Signal Preprocessing
- Independent Component Analysis (ICA)
- Artifact Subspace Reconstruction (ASR)
- Time-Frequency Analysis

## Key Takeaways

This project provided practical experience in biomedical signal processing and EEG data analysis.

The main learning outcomes include:

1. Understanding EEG frequency characteristics.
2. Processing and cleaning multichannel EEG recordings.
3. Applying STFT to non-stationary biological signals.
4. Extracting frequency-specific neural activity.
5. Analyzing stimulus-related changes in EEG power.
6. Comparing temporal neural responses across cognitive conditions.
7. Understanding the potential role of neural oscillations in cognitive impairment.

## Limitations

The study has several limitations:

- Only one subject per cognitive group was analyzed.
- The findings are exploratory rather than statistically validated.
- EEG preprocessing choices may influence the extracted power characteristics.
- Variability between individual subjects cannot be separated from group-level effects.
- Direct phase-amplitude coupling was not calculated in this phase.

Larger datasets and statistical validation would be required before drawing clinical conclusions.

## Future Work

Potential extensions include:

- Direct theta-gamma phase-amplitude coupling analysis.
- Modulation Index estimation.
- Baseline-normalized time-frequency analysis.
- Statistical comparisons across larger subject groups.
- Investigation of additional EEG channels and brain regions.
- Evaluation of EEG-derived features for cognitive impairment classification.



## Author

**Parnia Parhizi**

Electrical Engineering Undergraduate  
Sharif University of Technology

Signals and Systems Project — Spring 2025
