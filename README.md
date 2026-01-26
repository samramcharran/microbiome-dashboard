# Microbiome Dataset Discovery Dashboard

A Streamlit application for searching, scoring, and analyzing microbiome datasets from NCBI Sequence Read Archive (SRA).

## Features

- **SRA Search**: Query NCBI SRA for microbiome datasets using flexible search terms
- **Quality Scoring**: Automatic scoring based on:
  - Sequencing depth (read count)
  - Average read length
  - Metadata completeness
  - Sequencing platform
- **Interactive Data Table**: Filter and sort results by quality grade and platform
- **Quality Charts**: Visualize quality distributions, platform breakdown, and score components
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

The dashboard will open in your browser at `http://localhost:8501`.

### Search Examples

- `microbiome[All Fields] AND human[Organism]` - Human microbiome studies
- `16S[All Fields] AND gut[All Fields]` - 16S gut microbiome studies
- `metagenome[All Fields] AND soil[All Fields]` - Soil metagenomics
- `amplicon[All Fields] AND mouse[Organism]` - Mouse amplicon studies

## Quality Scoring System

Datasets are scored on a 100-point scale:

| Component | Max Points | Criteria |
|-----------|------------|----------|
| Sequencing Depth | 30 | Based on total read count |
| Read Length | 20 | Average read length |
| Metadata | 30 | Number of sample attributes |
| Platform | 20 | Sequencing technology used |

Grades are assigned as:
- **A**: 90-100 points
- **B**: 80-89 points
- **C**: 70-79 points
- **D**: 60-69 points
- **F**: Below 60 points

## Project Structure

```
microbiome-dashboard/
├── app.py              # Main Streamlit application
├── requirements.txt    # Python dependencies
└── README.md          # Documentation
```

## Technical Details

- Uses NCBI E-utilities API for SRA queries
- Parses SRA XML responses for metadata extraction
- Built with Streamlit for interactive web interface
- Plotly for responsive visualizations
- Pandas for data manipulation

## License

MIT License
