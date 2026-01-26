# Microbiome Dataset Discovery Dashboard

A Streamlit application for searching, scoring, and analyzing microbiome datasets from NCBI Sequence Read Archive (SRA). Optimized for both short-read (Illumina) and long-read (Oxford Nanopore) sequencing technologies.

## Features

- **SRA Search**: Query NCBI SRA for microbiome datasets using flexible search terms
- **Technology-Specific Scoring**: Separate scoring algorithms for:
  - **Short-read (Illumina)**: Optimized for read count and 250bp+ reads
  - **Long-read (Nanopore)**: Optimized for throughput (Gb) and read length (10kb+)
- **Clinical Metadata Detection**: Identifies datasets with patient/clinical metadata relevant for large cohort studies
- **Fecal Sample Filtering**: Filter for gut microbiome-relevant stool/fecal samples
- **Interactive Data Table**: Filter and sort results by quality grade, sequencing type, and sample type
- **Quality Charts**: Visualize quality distributions, platform breakdown, and score components
- **Nanopore Analysis Tab**: Dedicated view for long-read data with read length vs throughput scatter plots
- **CSV Export**: Download filtered results for further analysis

## Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/microbiome-dashboard.git
cd microbiome-dashboard

# Create virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

### Search Examples

**Nanopore Fecal Studies:**
```
fecal[All Fields] AND Oxford Nanopore[Platform]
```

**Long-Read Gut Microbiome:**
```
gut microbiome[All Fields] AND (nanopore[All Fields] OR "long read"[All Fields])
```

**Clinical Stool Studies:**
```
stool[All Fields] AND clinical[All Fields] AND microbiome[All Fields]
```

**Human Fecal 16S:**
```
16S[All Fields] AND fecal[All Fields] AND human[Organism]
```

## Quality Scoring System

Datasets are scored on a 100+ point scale with technology-specific criteria:

### Short-Read (Illumina) Scoring

| Component | Max Points | Criteria |
|-----------|------------|----------|
| Sequencing Depth | 30 | Based on total read count |
| Read Length | 20 | 250bp+ scores highest |
| Metadata | 30 | Number of sample attributes |
| Platform | 20 | Illumina scores highest |
| Clinical Fields | 10 | Patient/clinical metadata |
| Sample Relevance | 5 | Fecal sample bonus |

### Long-Read (Nanopore) Scoring

| Component | Max Points | Criteria |
|-----------|------------|----------|
| Throughput | 30 | Based on total gigabases (10+ Gb optimal) |
| Read Length | 25 | 10kb+ scores highest (strain-level resolution) |
| Metadata | 25 | Number of sample attributes |
| Platform | 20 | PromethION > MinION/GridION > Flongle |
| Clinical Fields | 10 | Patient/clinical metadata |
| Sample Relevance | 5 | Fecal sample bonus |

### Quality Grades

- **A**: 90+ points
- **B**: 80-89 points
- **C**: 70-79 points
- **D**: 60-69 points
- **F**: Below 60 points

### Clinical Metadata Fields Tracked

Subject ID, patient ID, collection date, timepoint, visit, age, sex, disease status, diagnosis, treatment, medication, BMI, diet, antibiotics, geographic location

## Project Structure

```
microbiome-dashboard/
├── app.py              # Main Streamlit application
├── requirements.txt    # Python dependencies
└── README.md          # Documentation
```

## Data Sources & References

All data is retrieved from official NCBI repositories. Every accession number is real and verifiable:

- **[NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)** - Primary source for sequencing datasets
- **[NCBI E-utilities API](https://www.ncbi.nlm.nih.gov/books/NBK25501/)** - Programmatic data retrieval
- **[NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/)** - Study and project metadata
- **[PubMed](https://pubmed.ncbi.nlm.nih.gov/)** - Linked scientific publications

### Accession Types

| Accession | Example | Description |
|-----------|---------|-------------|
| SRR | SRR36753138 | Sequencing run (unique per dataset) |
| SRX | SRX31751964 | Experiment accession |
| PRJNA | PRJNA1399953 | BioProject (groups related samples) |
| PMID | 12345678 | PubMed publication ID |

## Technical Details

- Uses NCBI E-utilities API for SRA queries
- Parses SRA XML responses for metadata extraction
- Built with Streamlit for interactive web interface
- Plotly for responsive visualizations
- Pandas for data manipulation

## Use Cases

- **Benchmarking**: Compare your sequencing data quality against public datasets
- **Reference Selection**: Find high-quality datasets for method validation
- **Clinical Study Design**: Identify comparable large-cohort studies
- **Technology Evaluation**: Compare Nanopore vs Illumina dataset characteristics

## License

MIT License
