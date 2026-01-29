# Microbiome Dataset Discovery Dashboard

**Live Demo:** https://microbiome-dashboard-samramcharran.streamlit.app

A Streamlit application for searching, scoring, and analyzing microbiome datasets from NCBI Sequence Read Archive (SRA). Optimized for both short-read (Illumina) and long-read (Oxford Nanopore) sequencing technologies. Built for genome banking, metadata curation, and strategic data prioritization workflows.

## Features

### Core Discovery Features
- **SRA Search**: Query NCBI SRA for microbiome datasets using flexible search terms
- **Technology-Specific Scoring**: Separate scoring algorithms for short-read (Illumina) and long-read (Nanopore)
- **Clinical Metadata Detection**: Identifies datasets with patient/clinical metadata relevant for large cohort studies
- **Fecal Sample Filtering**: Filter for gut microbiome-relevant stool/fecal samples
- **Disease Classification**: Auto-categorizes datasets by disease area (Gut-Brain Axis, Pain, GI Disorders, etc.)

### Mission Control Dashboard (Leadership View)
- **Strategic Metrics**: Genomes identified, banked, pending review, gut-brain relevant, total data volume
- **Disease Burden Priority Matrix**: Visualizes strategic value of disease categories with weighted scoring
- **Team Readiness Indicators**: Status for Wet Lab, Computational, and Clinical Analysis teams
- **Infrastructure Planning**: Storage, download time, and processing estimates
- **Strategic Recommendations**: Auto-generated action items based on current data

### Dataset Browser Enhancements
- **Favorites/Watchlist**: Star datasets for tracking, filter by favorites, export watchlist
- **Dataset Comparison**: Side-by-side comparison of up to 3 datasets
- **Date Filtering**: Filter datasets by collection date range
- **Interactive Data Table**: Sortable, filterable table with direct NCBI links

### Visualization Features
- **Temporal Trend Chart**: Bar chart of datasets by collection year
- **Platform Breakdown**: Pie chart of sequencing technologies and data volume by status
- **Discovery Funnel**: Visual pipeline from discovered to high-quality to public-ready
- **Disease Distribution Charts**: Bar and pie charts for disease categories

### Search URL Sharing
- **Shareable URLs**: Search parameters encoded in URL (q, role, preset, max, scoring)
- **Auto-Search**: URLs with parameters automatically execute search on load
- **Copy Search Link**: Easy sharing of search configurations

### Export Options
- **CSV Downloads**: Full results, banked only, or watchlist
- **Accession Lists**: Plain text SRR IDs for bulk download
- **JSON Export**: Structured data for pipeline integration
- **Storage Estimates**: S3 cost calculations for banked data

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

### NCBI API Compliance

This application **fully complies with NCBI E-utilities guidelines** out of the box:

- ✅ **Tool identification**: All API requests include `tool=microbiome-dashboard`
- ✅ **Rate limiting**: Automatic 0.34-second delay between requests (stays under 3 requests/second limit)
- ✅ **No data hoarding**: Real-time queries only, no caching or mirroring of NCBI data

**Optional enhancements** (not required for compliance):

```bash
# Developer contact email (for NCBI service notifications)
export NCBI_EMAIL="your-email@example.com"

# API key for higher rate limits - 10/sec instead of 3/sec
# Get one at https://www.ncbi.nlm.nih.gov/account/settings/
export NCBI_API_KEY="your-api-key"
```

Or create a `.env` file in the project root (not committed to git).

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

### URL Parameters for Search Sharing

Share searches by appending parameters to your dashboard URL:

```
?q=fecal[All Fields] AND microbiome[All Fields]&role=Leadership&preset=Custom Query&max=50&scoring=auto
```

Parameters:
- `q`: Search query (NCBI Entrez syntax)
- `role`: Role preset (Leadership, Wet Lab, Computational, BD / Partnerships, Disease Focus)
- `preset`: Preset query name or "Custom Query"
- `max`: Maximum results (10-100)
- `scoring`: Scoring mode (auto, nanopore, illumina)

## Tab Structure

| Tab | Description |
|-----|-------------|
| **Overview** | Discovery funnel, quick insights, temporal trend chart, platform breakdown |
| **Mission Control** | Leadership dashboard with strategic metrics, disease burden matrix, team readiness, and recommendations |
| **Dataset Browser** | Filterable table with favorites, comparison, date filtering |
| **Disease Cohorts** | Disease category distribution, gut-brain relevance, cohort prioritization |
| **Data Quality** | Vault status, access type, sequencing technology, sample type distributions |
| **Export** | Download options for CSV, JSON, accession lists, and watchlist |

## How It Works

### Two Views

The application has two main modes, accessible via a toggle in the sidebar:

| View | Purpose | Who it's for |
|------|---------|--------------|
| **Dataset Discovery** | Search and browse microbiome datasets | Scientists, computational biologists, researchers |
| **Analytics Dashboard** | Strategic metrics and operations overview | Leadership, executives, strategic planning |

### Dataset Discovery (Default)

The main landing experience for discovering microbiome data:

| Search Category | Options |
|-----------------|---------|
| **By Condition** | Mental Health, Pain Conditions, Digestive Health, Metabolic Health |
| **By Data Type** | Long-Read Sequencing, Shotgun Metagenomics, 16S Amplicon |
| **By Study Design** | Large Cohorts, Clinical Trials |
| **Browse** | All Available Data, Custom Search |

**Tabs:** Browse Datasets, Overview, Disease Categories, Quality Details, Export

### Analytics Dashboard (Leadership)

Strategic operations dashboard with focused views:

| Focus Area | What It Shows |
|------------|---------------|
| **Banking Progress** | Live status of identified vs banked genomes, pipeline funnel |
| **Disease Prioritization** | Strategic disease burden matrix, category breakdown |
| **Data Quality** | Quality score distribution, scoring factors, platform comparison |
| **Access Status** | Public vs private data, partnership opportunities |
| **Infrastructure** | Storage requirements, download estimates, team readiness |

**Sections:** Mission Control metrics at top, detailed views in tabs below

## Quality Scoring System

Datasets are scored on a 100+ point scale with technology-specific criteria:

### Short-Read (Illumina) Scoring

| Component | Max Points | Criteria |
|-----------|------------|----------|
| Sequencing Depth | 25 | Based on total read count |
| Read Length | 20 | 250bp+ scores highest |
| Metadata | 25 | Number of harmonized sample attributes |
| Clinical Fields | 15 | Patient/clinical metadata |
| Sample Relevance | 10 | Fecal sample bonus |
| Publication | 5 | Linked PubMed ID |

### Long-Read (Nanopore) Scoring

| Component | Max Points | Criteria |
|-----------|------------|----------|
| Throughput | 25 | Based on total gigabases (10+ Gb optimal) |
| Read Length | 25 | 10kb+ scores highest (strain-level resolution) |
| Metadata | 25 | Number of harmonized sample attributes |
| Clinical Fields | 15 | Patient/clinical metadata |
| Sample Relevance | 10 | Fecal sample bonus |
| Publication | 5 | Linked PubMed ID |

### Vault Status Thresholds

- **Banked**: 70+ points - High-quality, ready for genome vault
- **Pending Review**: 55-69 points - Requires manual review
- **Not Suitable**: <55 points - Does not meet quality thresholds

### Disease Burden Weights

Strategic prioritization weights for Mission Control:

| Category | Weight |
|----------|--------|
| Gut-Brain Axis | 10 |
| Pain Conditions | 9 |
| GI Disorders | 8 |
| Metabolic | 7 |
| Immune/Inflammatory | 6 |
| Infectious | 5 |
| Healthy/Control | 3 |
| Unclassified | 1 |

## Project Structure

```
microbiome-dashboard/
├── app.py              # Main Streamlit application
├── requirements.txt    # Python dependencies
└── README.md          # Documentation
```

## Dependencies

- streamlit>=1.28.0
- pandas>=2.0.0
- plotly>=5.18.0
- requests

## Disclaimer

**This application is for educational and research exploration purposes only.** It is not intended to replace, replicate, or serve as a substitute for NCBI's official services. Users should always verify data directly at [NCBI](https://www.ncbi.nlm.nih.gov/) for research and publication purposes.

## Data Sources & References

All data is retrieved from official NCBI repositories via their public API. Every accession number is real and verifiable:

- **[NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)** - Primary source for sequencing datasets
- **[NCBI E-utilities API](https://www.ncbi.nlm.nih.gov/books/NBK25497/)** - Programmatic data retrieval
- **[NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/)** - Study and project metadata
- **[PubMed](https://pubmed.ncbi.nlm.nih.gov/)** - Linked scientific publications

### NCBI Data Usage

NCBI databases are produced by a U.S. government agency and are in the public domain. This application complies with [NCBI's E-utilities usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/):
- Requests include tool identification and respect rate limits
- No bulk downloading or mirroring of NCBI data
- Data is retrieved in real-time and not cached or redistributed

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
- Plotly for responsive visualizations including scatter_geo maps
- Pandas for data manipulation
- Session state for favorites, comparison selections, and search results

## Use Cases

- **Genome Banking**: Identify high-quality datasets for vault ingestion
- **Strategic Planning**: Use Mission Control for leadership decision-making
- **Wet Lab Coordination**: Find long-read datasets for strain isolation
- **Clinical Analysis**: Prioritize gut-brain relevant datasets
- **Partnership Outreach**: Identify restricted-access datasets for collaboration
- **Benchmarking**: Compare your sequencing data quality against public datasets
- **Technology Evaluation**: Compare Nanopore vs Illumina dataset characteristics

## License

MIT License

**Note:** This license applies to the application code only. Data retrieved from NCBI is in the public domain as a work of the U.S. Government. NCBI requests that users cite the original data sources in publications. See [NCBI's policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) for details.
