<a name="readme-top"></a>

<h3 align="center">DNA Methylation Analysis Pipeline</h3>

  <p align="center">
   Python-based analysis pipeline to process raw Next-Generation Sequencing data and produce a comprehensive analysis of the desired locus.
    <br />
    <a href="https://github.com/pablo-romeroclavijo/DNA_Methylation_Analysis_Pipeline/tree/main#readme"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    ·
    <a href="https://github.com/pablo-romeroclavijo/Pokedex_front_end/issues">Report Bug</a>
    ·
    <a href="https://github.com/pablo-romeroclavijo/Pokedex_front_end/issues">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
        <li><a href="#flow-chart">Flow-Chart</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Pre-requisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>

  </ol>
</details>

<!-- ABOUT THE PROJECT -->

## About The Project

<!-- <img src="https://github.com/LaFosseAcademy/debug-assignment-2-pablo-romeroclavijo/assets/136720541/9f8670c0-6f18-4135-a012-5e84fc23c5cf" alt="Demo" width="150" height="150"> -->

This analysis pipeline is design to process DNA methylation results originated from Illumina Next Generation Sequencing. Approximately 20.000 sequences at the time.

The programme is divided in 3 sections:

    1.Decoder: demultiplexes the data into subsamples based on barcoding of the DNA sequencing.  Data is stores in a new directory 1_Decoded_data. Runtime: 2min approx.

    2. Analyser: for every subsample, it will filter and compare the individual sequence with a provided reference. It will calculate several biologically relevant parameter and correct for technique bias. Data is stored on individual data files 2_Methylation results.

    3. Plotting: for this example very basic plots have been included. Including QC plots. Data stored in new Directory 3_Plot.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With

![npm package](https://img.shields.io/badge/Python-v3.6-brightgreen.svg)
![npm package](https://img.shields.io/badge/Pandas-v2.1-brightgreen.svg)
![npm package](https://img.shields.io/badge/Numpy-v1.26.1-brightgreen.svg)
![npm package](https://img.shields.io/badge/MatPlotLib-v3.7.1-brightgreen.svg)
![npm package](https://img.shields.io/badge/BioPython-v1.81-brightgreen.svg)

### Flow Chart

<img src="https://github.com/pablo-romeroclavijo/DNA_Methylation_Analysis_Pipeline/blob/main/WorkFLow.jpg?raw=true" alt = "flowchart" width=1000px>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->

## Getting Started

### Pre-requisites

Download lastest Python and Pip versions (<a href="https://www.python.org/downloads/macos/">Python > Downloads > MacOS</a>)

  ```sh
    python -m ensurepip --upgrade
  ```

### Installation

1. Clone the repo

   ```sh
      git clone https://github.com/pablo-romeroclavijo/DNA_Methylation_Analysis_Pipeline.git
   ```

2. Install dependencies
   ```sh
      pip install requirements.txt
   ```
3. Open the Decoder+Quantifier.py and decide what parts of the pipeline to run. (Comment out line 556 or 563)

   **Code Block**

```py
    554          ## RUN Decoder
    555          print(' ---------------------------  \n DECODER ')
    556          decode = Decoder(working_directory)
    557          decode.run()
    558          print('\n ---------------------------\n')
    559
    560          ## Run Quantification and plotting
    561          print('\n ---------------------------  \n QUANTIFIER')
    562
    563          quantify = Quantifier(working_directory)
    564          quantify.run()
```

4. Run the entire application

   ```sh
      python Decoder+Quantifier.py
   ```

5. Select data directory(default = <"current-directory">)

<!-- USAGE EXAMPLES -->

## Required files and directories
<ol>
<li>0_Raw_data:</li>
  <ul>
    <li> Fastq files: can be modified in the SeqIO parameters in case of using FASTA</li>
  </ul>
<li>0_Reference</li>
  <ul>
    <li>Ref_uncon.txt ---- Fasta sequence of the unconverted amplicon</li>
    <li>Barcodes.txt ---- file required for Decoder</li>
  </ul>
</ol>

<!-- CONTACT -->

## Contact

Pablo Romero - [@Github](https://github.com/pablo-romeroclavijo) // [@LinkedIn](https://www.linkedin.com/in/pablo-romeroclavijo/)

Project Link: [https://github.com/pablo-romeroclavijo/DNA_Methylation_Analysis_Pipeline](https://github.com/pablo-romeroclavijo/DNA_Methylation_Analysis_Pipeline/tree/main#readme)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
